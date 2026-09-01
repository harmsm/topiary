import pytest
import topiary
from topiary.quality import remove_redundancy
from topiary.quality.redundancy import _get_quality_scores
from topiary.quality.redundancy import _construct_args
from topiary.quality.redundancy import _compare_seqs
from topiary.quality.redundancy import _redundancy_thread_function
from topiary.quality.redundancy import _EXPECTED_COLUMNS
from topiary.quality.redundancy import find_redundancy_cutoff

import numpy as np
import pandas as pd
import multiprocessing

@pytest.mark.smoke
def test__get_quality_scores(test_dataframes):

    # Get copy of the dataframe -- we're going to hack it
    df = test_dataframes["good-df"].copy()
    species_in_df = list(df.species)

    # Needs to be defined.
    df["diff_from_median"] = 0

    # Make sure that the key species encoding works as expected
    # No key_species passed in or key_species not in dataframe
    assert _get_quality_scores(df.loc[0,:])[1] == 1
    assert _get_quality_scores(df.loc[0,:],{"Not a species":None})[1] == 1

    # Key species
    assert _get_quality_scores(df.loc[0,:],{species_in_df[0]:None})[1] == 0

    # Make sure quality assignment doing what we think
    for i in df.index:
        row = df.loc[i,:]
        scores = _get_quality_scores(row)

        values_from_df = list(row[_EXPECTED_COLUMNS])
        values_from_df.append(1/len(row.sequence))
        values_from_df = np.array(values_from_df,dtype=float)

        assert np.array_equal(scores[2:],values_from_df)

    # Test always_keep
    # No always_keep in datafframe
    assert _get_quality_scores(df.loc[0,:])[0] == 1

    df["always_keep"] = False # always keep False
    assert _get_quality_scores(df.loc[0,:])[0] == 1

    df["always_keep"] = True # always keep True
    assert _get_quality_scores(df.loc[0,:])[0] == 0

    # Test ValueError for target_length without pct_length_cutoff
    with pytest.raises(ValueError):
        _get_quality_scores(df.loc[0,:], target_length=100)

def test__construct_args():

    sequence_array = np.array(["STARE" for _ in range(4)])
    quality_array = np.array([np.zeros(3,dtype=int) for _ in range(4)])
    keep_array = np.ones(4,dtype=bool)
    cutoff = 0.9
    discard_key = True

    kwargs_list, num_threads = _construct_args(sequence_array=sequence_array,
                                               quality_array=quality_array,
                                               keep_array=keep_array,
                                               cutoff=cutoff,
                                               discard_key=discard_key,
                                               num_threads=1)

    assert num_threads == 1
    assert np.array_equal(kwargs_list[0]["i_block"],(0,4))
    assert np.array_equal(kwargs_list[0]["j_block"],(0,4))

    out = []
    for a in kwargs_list:
        i_block = a["i_block"]
        j_block = a["j_block"]
        for i in range(i_block[0],i_block[1]):
            for j in range(j_block[0],j_block[1]):
                if i >= j:
                    continue
                out.append((i,j))

    out.sort()
    assert np.array_equal(out,((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)))


    # Not worth chopping up problem for this small of an array -- should set
    # number of threads to 1.
    kwargs_list, num_threads = _construct_args(sequence_array=sequence_array,
                                               quality_array=quality_array,
                                               keep_array=keep_array,
                                               cutoff=cutoff,
                                               discard_key=discard_key,
                                               num_threads=2)

    assert num_threads == 1
    assert np.array_equal(kwargs_list[0]["i_block"],(0,4))
    assert np.array_equal(kwargs_list[0]["j_block"],(0,4))

@pytest.mark.smoke
def test__compare_seqs(test_dataframes):

    A_seq = "TEST"
    B_seq = "TAST"

    # Identical quals
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)

    # Neither are always keep sequences
    A_qual[0] = 1
    B_qual[0] = 1

    # Neither are key sequences
    A_qual[1] = 1
    B_qual[1] = 1

    # Keep both; below cutoff
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.9)
    assert a1 is True
    assert a2 is True

    # Keep A arbitrarily
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is False

    # Now make A_qual score worse than B, so keep B
    A_qual[2] = 1
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # Not set up qual scores so neither are key_species, B has earlier better
    # score than A
    A_qual = np.ones(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.ones(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[-1] = 0
    B_qual[-2] = 0

    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # both key species, A worse than B. No always_keep
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[0] = 1
    B_qual[0] = 1
    A_qual[2] = 1

    # implicit discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is True

    # Explicit discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5,discard_key=False)
    assert a1 is True
    assert a2 is True

    # Check discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5,discard_key=True)
    assert a1 is False
    assert a2 is True

    # both always keep, but  beter. Should keep both
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[2] = 1
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is True

def test__redundancy_thread_function():
    # Setup sequences and quality
    sequence_array = np.array(["AAAA", "AAAA", "BBBB"])
    quality_array = np.array([[1,1,0,0,0,0,0,0,0], [1,1,0,0,0,0,0,0,0], [1,1,0,0,0,0,0,0,0]])
    keep_array = multiprocessing.Array('i', [1, 1, 1])
    lock = multiprocessing.Lock()
    
    # Test redundancy within block (0, 3) and (0, 3)
    _redundancy_thread_function((0, 3), (0, 3), sequence_array, quality_array, keep_array, 0.5, False, lock)
    
    # AAAA vs AAAA should drop one. 
    # BBBB is unique.
    assert keep_array[0] == 1
    assert keep_array[1] == 0
    assert keep_array[2] == 1


@pytest.mark.smoke
def test_remove_redundancy(test_dataframes):

    df = test_dataframes["good-df"].copy()

    # -------------------------------------------------------------------------
    # Test argument parsing

    bad_df = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_df:
        with pytest.raises(ValueError):
            remove_redundancy(df=b)

    remove_redundancy(df=df)

    bad_cutoff = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_cutoff:
        with pytest.raises(ValueError):
            remove_redundancy(df=df,cutoff=b)

    good_cutoff = [0,0.5,1]
    for g in good_cutoff:
        remove_redundancy(df=df,cutoff=g)

    bad_silent = [None,"test",int,float,{"test":1}]
    for b in bad_silent:
        with pytest.raises(ValueError):
            remove_redundancy(df=df,silent=b)

    good_silent = [True,False,0,1]
    for g in good_silent:
        remove_redundancy(df=df,silent=g)


    # -------------------------------------------------------------------------
    # Make sure dropping is happening a sane way that depends on cutoff and
    # key_species.

    df = test_dataframes["good-df"].copy()
    species_in_df = list(df.species)

    # sequences in this dataframe are between 0.9125 and 0.98125 identical.
    out_df = remove_redundancy(df=df,cutoff=0.99)
    assert np.sum(out_df.keep) == np.sum(df.keep)

    # Cut some
    out_df = remove_redundancy(df=df,cutoff=0.96)
    assert np.sum(out_df.keep) < np.sum(df.keep)

    # Cut basically all -- only one shoudl survive
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert np.sum(out_df.keep) == 1

    # All key species -- all should survive
    df["key_species"] = True
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert np.sum(out_df.keep) == np.sum(df.keep)

    # One isn't keep -- make sure it's dropped
    df.loc[0,"key_species"] = False
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert np.sum(out_df.keep) == 4
    assert out_df.loc[out_df["species"] == species_in_df[0],:].iloc[0].keep == False

    # Test mismatching column types
    df_bad_col = df.copy()
    df_bad_col["partial"] = ["not a number"] * len(df)
    with pytest.raises(ValueError):
        remove_redundancy(df_bad_col)

    # -------------------------------------------------------------------------
    # Make sure it takes row with higher quality

    df = test_dataframes["good-df"].copy()
    df = df.iloc[2:4]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == False
    assert out_df.keep.iloc[1] == True

    df = test_dataframes["good-df"].copy()
    df = df.iloc[1:3]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == True
    assert out_df.keep.iloc[1] == False

    # Make sure length bit is being processed properly. Make first sequence short
    # so it gets dropped
    df = test_dataframes["good-df_only-required-columns"].copy()
    df.loc[0,"sequence"] = "MLPFLFFS"
    df.loc[:,"length"] = [len(s) for s in df.loc[:,"sequence"]]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == False

    # -------------------------------------------------------------------------
    # Check input dataframe without quality information

    # Should work fine but cut nothing
    df = test_dataframes["good-df_only-required-columns"].copy()
    out_df = remove_redundancy(df=df,cutoff=0.99)
    assert np.sum(out_df.keep) == len(df.sequence)

    # Cut some
    out_df = remove_redundancy(df=df,cutoff=0.96)
    assert np.sum(out_df.keep) <= len(df.sequence)

    # Cut basically all -- only one should survive
    out_df = remove_redundancy(df=df,cutoff=0.2)
    assert np.sum(out_df.keep) == 1

@pytest.mark.smoke
def test_find_redundancy_cutoff(test_dataframes):
    df = test_dataframes["good-df"].copy()
    
    # Test target reached immediately (max_cutoff)
    # len(df) is small, so we'll likely get a quick result
    cutoff = find_redundancy_cutoff(df, target_seq_number=10, max_cutoff=0.99)
    assert 0.25 <= cutoff <= 0.99
    
    # Test target reached immediately (min_cutoff)
    cutoff = find_redundancy_cutoff(df, target_seq_number=1, min_cutoff=0.25)
    assert 0.25 <= cutoff <= 0.99
    
    # Test optimization loop
    # We'll use a larger target than 1 to force some iterations
    # But good-df is only 5 sequences.
    # Let's mock a larger df if needed, or just test current logic.
    cutoff = find_redundancy_cutoff(df, target_seq_number=3)
    assert 0.25 <= cutoff <= 0.99
    
    # Test min > max raises ValueError
    with pytest.raises(ValueError):
        find_redundancy_cutoff(df, target_seq_number=3, min_cutoff=0.9, max_cutoff=0.8)

    # Test min == max returns immediately
    assert find_redundancy_cutoff(df, target_seq_number=3, min_cutoff=0.5, max_cutoff=0.5) == 0.5

def _graded_similarity_df(n=10):
    """
    Sequences with graded divergence, so different identity cutoffs give
    genuinely different keep counts. Needed to make find_redundancy_cutoff
    actually search rather than hitting one of its early returns.
    """

    base = "MASTPDLLKWAQRSTVYNEGHIKLMNPQRSTV"

    seqs = []
    for i in range(n):
        s = list(base)
        for j in range(i):
            s[j] = "A"
        seqs.append("".join(s))

    return pd.DataFrame({
        "name":[f"s{i}" for i in range(n)],
        "species":[f"Species {i}" for i in range(n)],
        "sequence":seqs,
        "keep":[True]*n,
        "uid":[f"uid{i:07d}" for i in range(n)],
    })


@pytest.mark.smoke
def test_find_redundancy_cutoff_searches_between_bounds():
    """
    The interesting path: the target lies strictly between what max_cutoff and
    min_cutoff produce, so the function has to run its binary search.

    sample_fx=1.0 so the whole dataframe is used and the result is
    deterministic (the default samples a random half).
    """

    df = _graded_similarity_df()

    # cutoff 0.99 keeps 9, cutoff 0.25 keeps 1 -- so a target of 3 is reachable
    # only by searching.
    cutoff = find_redundancy_cutoff(df,
                                    target_seq_number=3,
                                    sample_fx=1.0,
                                    max_cutoff=0.99,
                                    min_cutoff=0.25,
                                    num_threads=1)

    assert 0.25 <= cutoff <= 0.99

    # And the cutoff it found actually gets near the target
    out = remove_redundancy(df.copy(),cutoff=cutoff,silent=True,num_threads=1)
    assert abs(np.sum(out.keep) - 3) <= 2


@pytest.mark.smoke
def test_find_redundancy_cutoff_respects_max_iterations():
    """
    With a tiny iteration budget the search stops early and returns the best
    cutoff seen so far rather than looping forever.
    """

    df = _graded_similarity_df()

    cutoff = find_redundancy_cutoff(df,
                                    target_seq_number=3,
                                    sample_fx=1.0,
                                    max_cutoff=0.99,
                                    min_cutoff=0.25,
                                    max_iterations=1,
                                    num_threads=1)

    assert 0.25 <= cutoff <= 0.99


@pytest.mark.smoke
def test_find_redundancy_cutoff_exact_hit_returns_immediately():
    """
    If a probed cutoff lands exactly on the target, that cutoff is returned
    without further searching.
    """

    df = _graded_similarity_df()

    # cutoff 0.99 keeps exactly 9, so asking for 9 returns max_cutoff directly
    cutoff = find_redundancy_cutoff(df,
                                    target_seq_number=9,
                                    sample_fx=1.0,
                                    max_cutoff=0.99,
                                    min_cutoff=0.25,
                                    num_threads=1)

    assert cutoff == 0.99
