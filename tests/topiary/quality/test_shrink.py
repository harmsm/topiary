
import pytest

import topiary
from topiary.quality.shrink import shrink_in_species
from topiary.quality.shrink import shrink_redundant
from topiary.quality.shrink import shrink_aligners
from topiary.quality.shrink import shrink_dataset


import numpy as np
import pandas as pd

import os

def test_shrink_in_species():
    """
    Drops near-identical sequences *within* a species, ignoring paralog calls.
    No species tree involved, so this is hermetic.
    """

    # Two species. Within Homo sapiens, seq A and A' are identical and should
    # collapse to one. B is different and survives. Mus musculus has a single
    # sequence identical to a human one -- but it is a different species, so it
    # must NOT be merged away.
    df = pd.DataFrame({
        "name":["hA","hA_dup","hB","mA"],
        "species":["Homo sapiens","Homo sapiens","Homo sapiens","Mus musculus"],
        "sequence":["MASTPDLLKW","MASTPDLLKW","QQQQQYYYYY","MASTPDLLKW"],
        "keep":[True,True,True,True],
        "uid":["aaaaaaaaaa","bbbbbbbbbb","cccccccccc","dddddddddd"],
    })

    out = shrink_in_species(df,redundancy_cutoff=0.98)

    kept = set(out.loc[out.keep,"uid"])

    # Exactly one of the two identical human sequences survives
    assert len({"aaaaaaaaaa","bbbbbbbbbb"} & kept) == 1

    # The distinct human sequence survives
    assert "cccccccccc" in kept

    # The mouse sequence survives despite matching a human sequence -- merging
    # is within-species only
    assert "dddddddddd" in kept

    # Nothing is deleted, only flagged
    assert len(out) == len(df)


def test_shrink_in_species_respects_always_keep():

    df = pd.DataFrame({
        "name":["a","b"],
        "species":["Homo sapiens","Homo sapiens"],
        "sequence":["MASTPDLLKW","MASTPDLLKW"],
        "keep":[True,True],
        "uid":["aaaaaaaaaa","bbbbbbbbbb"],
        "always_keep":[False,True],
    })

    out = shrink_in_species(df,redundancy_cutoff=0.98)

    # The always_keep sequence survives even though it is redundant
    assert bool(out.loc[out["uid"] == "bbbbbbbbbb","keep"].iloc[0])


def test_shrink_in_species_validates_cutoff():

    df = pd.DataFrame({
        "name":["a"],
        "species":["Homo sapiens"],
        "sequence":["MASTPDLLKW"],
        "keep":[True],
        "uid":["aaaaaaaaaa"],
    })

    for bad in [-0.1,1.1,"not_a_float",None]:
        with pytest.raises(ValueError):
            shrink_in_species(df,redundancy_cutoff=bad)

@pytest.mark.network
def test_shrink_redundant(for_real_inference):

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])

    # Test with default parameters
    df2 = shrink_redundant(df)
    expected_keep = np.array(df.keep)
    expected_keep[1] = False
    assert np.array_equal(expected_keep,df2.keep)

    # check arg checking
    with pytest.raises(ValueError):
        shrink_redundant(df,paralog_column="not_a_column")

    with pytest.raises(ValueError):
        shrink_redundant(df,weighted_paralog_split="stupid")

    with pytest.raises(ValueError):
        shrink_redundant(df,merge_block_size=-1)

    with pytest.raises(ValueError):
        shrink_redundant(df,redundancy_cutoff=-1)

    # Construct test dataset with only two species
    test_df = df.loc[df.species.isin(["Xenopus laevis","Geotrypetes seraphini"]),:].copy()
    test_df.loc[:,"recip_paralog"] = ["LY96","LY96","LY96","LY86","LY86"]
    test_df.loc[:,"always_keep"] = np.zeros(len(test_df),dtype=bool)
    test_df.loc[:,"key_species"] = np.zeros(len(test_df),dtype=bool)

    # Should whack from LY86 and LY96
    out_df = shrink_redundant(test_df,redundancy_cutoff=0.5)
    assert np.array_equal(out_df.keep,np.array([False,True,True,True,False]))

    # Should whack from LY86 and LY96
    out_df = shrink_redundant(test_df,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,False,True,True,False]))

    # Should keep everything
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"always_keep"] = np.ones(len(test_df),dtype=bool)
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,True,True,True,True]))

    # Should keep only one now that everything is same paralog
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = ["LY96","LY96","LY96","LY96","LY96"]
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,False,False,True,False]))

    # Should keep everything because all are different paralogs
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = ["A","B","C","D","E"]
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,True,True,True,True]))

    # Make a sequence get not kept because it is very short
    test_df_2 = test_df.copy()
    test_df_2.loc[[9],"sequence"] = "MEVW"
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.5)
    assert np.array_equal(out_df.keep,np.array([False,False,True,True,False]))

    # All have identical sequence and paralog. Take first because it is the
    # best -- not low_quality
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = "LY96"
    test_df_2.loc[:,"sequence"] = test_df_2.loc[8,"sequence"]
    test_df_2.loc[:,"low_quality"] = True
    test_df_2.loc[8,"low_quality"] = False
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,False,False,False,False]))

    # All have identical sequence and paralog. Take second because it is not
    # low quality -- rest are.
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = "LY96"
    test_df_2.loc[:,"sequence"] = test_df_2.loc[8,"sequence"]
    test_df_2.loc[:,"low_quality"] = False
    test_df_2.loc[:,"partial"] = True
    test_df_2.loc[9,"partial"] = False
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,True,False,False,False]))

def _aligner_df():
    """
    Minimal dataframe with the columns shrink_aligners needs: two paralogs,
    four species, one flagged as a key species.
    """

    return pd.DataFrame({
        "name":["a","b","c","d"],
        "species":["Homo sapiens","Mus musculus","Gallus gallus","Danio rerio"],
        "sequence":["MASTPDLLKW","MASTPDLLKQ","QQQQQYYYYY","QQQQQYYYYW"],
        "keep":[True,True,True,True],
        "uid":["aaaaaaaaaa","bbbbbbbbbb","cccccccccc","dddddddddd"],
        "recip_paralog":["p1","p1","p2","p2"],
        "key_species":[True,False,False,False],
        "always_keep":[False,False,False,False],
    })


@pytest.mark.smoke
def test_shrink_aligners():
    """
    Trims a dataset down to a target size by alignment quality.

    species_tree_aware=False uses dummy merge blocks, so no species tree and no
    network -- but the alignment-quality scoring still shells out to muscle,
    which makes this smoke rather than unit. (The species-tree-aware path is
    covered by test_shrink_dataset, which is marked network.)
    """

    df = _aligner_df()

    out = shrink_aligners(df,target_seq_number=2,species_tree_aware=False)

    # Trimmed to the target, and nothing deleted -- only flagged
    assert sum(out.keep) == 2
    assert len(out) == len(df)


@pytest.mark.smoke
def test_shrink_aligners_at_full_target():
    """
    Asking for the whole dataset back still runs the alignment-quality step, so
    this one shells out to muscle -- hence smoke rather than unit.
    """

    out = shrink_aligners(_aligner_df(),
                          target_seq_number=4,
                          species_tree_aware=False)

    assert sum(out.keep) == 4


def test_shrink_aligners_honors_always_keep():

    df = _aligner_df()
    df["always_keep"] = [False,False,False,True]

    out = shrink_aligners(df,target_seq_number=1,species_tree_aware=False)

    # The always_keep sequence survives an aggressive target
    assert bool(out.loc[out["uid"] == "dddddddddd","keep"].iloc[0])


def test_shrink_aligners_validates_arguments():

    df = _aligner_df()

    for bad in ["not_a_bool",1.5,None]:
        with pytest.raises(ValueError):
            shrink_aligners(df,target_seq_number=2,species_tree_aware=bad)

@pytest.mark.network
@pytest.mark.skipif(os.name == "nt",reason="muscle cannot be installed via conda on windows")
def test_shrink_dataset(for_real_inference,tmpdir):

    # shrink_dataset calls align(), which writes topiary-tmp_* scratch files
    # into the working directory. Without this chdir they land in the repo root
    # -- and survive there if muscle raises before align() cleans up.
    os.chdir(tmpdir)

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])
    df.loc[df["recip_paralog"] == "unassigned","recip_paralog"] = "LY96"

    # Runs with default
    out = shrink_dataset(df)

    # no easy way to check that this *works* but at least this should not die
    # when other strings will
    out = shrink_dataset(df,paralog_column="nickname")

    # bad paralog_column values
    bad_paralog_column = ["not_a_column",None,14.2,str]
    for b in bad_paralog_column:
        print(b)
        with pytest.raises(ValueError):
            shrink_dataset(df,paralog_column=b)

    # This should not shrink at all. All sequences should make it fine, as
    # stringency is very low and we allow up to 0.12 seqs/per column. 160
    # columns in alignment, 19 sequences. 160*0.12 --> 19.2
    out = shrink_dataset(df,seqs_per_column=0.12,redundancy_cutoff=1.0,sparse_column_cutoff=1.0)
    assert np.sum(out.keep) == 19

    # Get only six sequence out, as we're only allowing 1/160 seqs per column.
    # There will be the four always_keep
    out = shrink_dataset(df,seqs_per_column=1/160)
    assert np.sum(out.keep) == 4

    tmp_df = df.copy()
    tmp_df.loc[:,"always_keep"] = False
    out = shrink_dataset(tmp_df,seqs_per_column=1/160)
    assert np.sum(out.keep) == 2

    # This will be only out.always_keep -- lowest we can go
    out = shrink_dataset(df,max_seq_number=1)
    assert np.sum(out.keep) == 4

    bad_max_seq_number = [-1,0,"test",int,14.2]
    for b in bad_max_seq_number:
        print(b)
        with pytest.raises(ValueError):
            shrink_dataset(df,max_seq_number=b)

    # IMPROVE TEST
    # THIS IS A PARTIAL TEST AT THE MOMENT. FOCUSING ON TESTING LOWER LEVEL
    # FUNCTIONS
