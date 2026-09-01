
import pytest

import topiary
from topiary._private import check

from topiary.io.seed import read_seed
from topiary.io.seed import df_from_seed

import numpy as np
import pandas as pd

import re


@pytest.mark.network
def test_read_seed(seed_dataframes,user_seed_dataframes):

    def _validate_output(out):

        df = out[0]
        key_species = out[1]
        paralog_patterns = out[2]
        species_aware = out[3]

        # check output dataframe
        assert len(df) == 8

        names = list(np.unique(df.name))
        names.sort()
        assert np.array_equal(names,("LY86","LY96"))

        species = list(np.unique(df.species))
        species.sort()
        assert np.array_equal(species,("Danio rerio",
                                       "Gallus gallus",
                                       "Homo sapiens",
                                       "Mus musculus"))

        assert np.array_equal(df.always_keep,np.ones(8,dtype=bool))

        # Make sure it's a good topiary dataframe
        check.check_topiary_dataframe(df)

        # Test key species
        assert np.array_equal(key_species,("Danio rerio",
                                           "Gallus gallus",
                                           "Homo sapiens",
                                           "Mus musculus"))

        # Test paralog patterns
        expected_patterns = {
            "LY96":"esop[\\ \\-_\\.]*1|ly[\\ \\-_\\.]*96|lymphocyte[\\ \\-_\\.]*antigen[\\ \\-_\\.]*96|myeloid[\\ \\-_\\.]*differentiation[\\ \\-_\\.]*protein[\\ \\-_\\.]*2",
            "LY86":"ly[\\ \\-_\\.]*86|lymphocyte[\\ \\-_\\.]*antigen[\\ \\-_\\.]*86|md[\\ \\-_\\.]*1|mmd[\\ \\-_\\.]*1|rp[\\ \\-_\\.]*105[\\ \\-_\\.]*associated[\\ \\-_\\.]*3"
        }
        for k in paralog_patterns:
            print(paralog_patterns[k].pattern)
            assert paralog_patterns[k].pattern == expected_patterns[k]

        assert species_aware is True

    # csv
    df_file = seed_dataframes["good-seed-df.csv"]
    out = read_seed(df_file)
    _validate_output(out)

    # Make sure we can read the dataframe as a dataframe object
    df_as_df = pd.read_csv(df_file)
    out = read_seed(df_as_df)
    _validate_output(out)

    # tsv
    df_file = seed_dataframes["good-seed-df.tsv"]
    out = read_seed(df_file)
    _validate_output(out)

    # xlsx
    df_file = seed_dataframes["good-seed-df.xlsx"]
    out = read_seed(df_file)
    _validate_output(out)

    # csv with .txt extension
    df_file = seed_dataframes["good-seed-df.txt"]
    out = read_seed(df_file)
    _validate_output(out)

    # bad df passes
    with pytest.raises(FileNotFoundError):
        read_seed("not-really-a-file")

    bad_df = [pd.DataFrame,None,1,[],{},1.5,str,int,float]
    for b in bad_df:
        print(f"passing bad_df {b}")
        with pytest.raises(ValueError):
            read_seed(b)

    # required column check
    good_df = pd.read_csv(seed_dataframes["good-seed-df.csv"])
    bad_df = good_df.drop(columns=["species"])
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # bad species
    bad_df = good_df.copy()
    bad_df.loc[:,"species"] = "not a species"
    with pytest.raises(ValueError):
        with pytest.warns():
            read_seed(bad_df)

    # species that is findable but not resolvable
    bad_df = good_df.copy()
    bad_df.loc[:,"species"] = "Bos indicus x Bos taurus"
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # species that is findable but not resolvable
    bad_df = seed_dataframes["duplicated-alias.xlsx"]
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # Read in a collection of user-generated dataframes, checking for expected
    # outputs
    
    is_species_aware = {'snase.xlsx':False,
                        'rnaseh.xlsx':False,
                        'chs.xlsx':True,
                        'iapp-cgrp.xlsx':True,
                        'gproteins.xlsx':False,
                        'ly86ly96.xlsx':True,
                        'zo1.xlsx': True,
                        's100a5-a6.xlsx':True}

    for s in user_seed_dataframes:

        print(f"Testing read of {s}")

        df, key_species, paralog_patterns, species_aware = topiary.io.read_seed(user_seed_dataframes[s])
        check_read = pd.read_excel(user_seed_dataframes[s])
        check_read.columns = [c.lower().strip() for c in check_read.columns]

        assert len(df) == len(check_read)

        for idx in df.index:
            assert df.loc[idx,"species"] == check_read.loc[idx,"species"].strip()
            assert df.loc[idx,"name"] == check_read.loc[idx,"name"].strip()

            this_seq = "".join(check_read.loc[idx,"sequence"].split())
            assert df.loc[idx,"sequence"] == this_seq

        assert np.array_equal(df.keep,np.ones(len(df),dtype=bool))
        assert np.array_equal(df.key_species,np.ones(len(df),dtype=bool))
        assert np.array_equal(df.always_keep,np.ones(len(df),dtype=bool))
        assert len(df.uid) == len(np.unique(df.uid))

        # Make sure the code is correctly identifying whether to treat this 
        # dataset as species aware
        assert species_aware is is_species_aware[s]



@pytest.fixture
def offline_seed(mocker,seed_dataframes):
    """
    df_from_seed reaches the outside world in four places: read_seed (Open Tree
    of Life), get_proteome_ids (NCBI), ott_to_mrca (OTL), and get_df_ott (OTL).
    Stub all four so the test exercises df_from_seed's own orchestration rather
    than the network.

    Returns the seed file path.
    """

    seed_file = seed_dataframes["good-seed-df.csv"]

    real_seed = pd.read_csv(seed_file)
    seed_df = real_seed.copy()
    seed_df["keep"] = True
    seed_df["always_keep"] = True
    seed_df["key_species"] = True
    seed_df["uid"] = [f"seed{i:06d}" for i in range(len(seed_df))]
    seed_df["ott"] = "ott9999"

    key_species = sorted(set(real_seed["species"]))
    paralog_patterns = {"LY96":[re.compile("lymphocyte antigen 96")],
                        "LY86":[re.compile("lymphocyte antigen 86")]}

    mocker.patch("topiary.io.read_seed",
                 return_value=(seed_df,key_species,paralog_patterns,True))
    mocker.patch("topiary.ncbi.entrez.get_proteome_ids",
                 return_value=("GCF_000001405.40",None))
    mocker.patch("topiary.opentree.ott_to_mrca",
                 return_value={"ott_name":"Vertebrata","ott_id":1234})

    def _fake_get_df_ott(df,**kwargs):
        df = df.copy()
        df["ott"] = "ott1234"
        df["orig_species"] = df["species"]
        df["phylo_context"] = "All life"
        return df

    mocker.patch("topiary.get_df_ott",side_effect=_fake_get_df_ott)

    # merge_and_annotate re-downloads sequences from Entrez for any hit it does
    # not fully trust. That call runs inside a multiprocessing pool, so on Linux
    # (fork) it inherits the network guard and dies with a confusing
    # PicklingError, while on macOS (spawn) it does not inherit the guard and
    # silently makes a real network call. Mock it so neither happens.
    mocker.patch("topiary.ncbi.blast.merge.get_sequences",
                 side_effect=lambda accessions,**kwargs: list(_sequences_for(accessions)))

    return seed_file


def _sequences_for(accessions):
    """A stand-in protein sequence per requested accession."""

    return ["MLPFLFFSTLFSSIFTEAQKQYWVCNSSDASISYTYCDKMQ" for _ in accessions]


def _real_blast_hits(ncbi_blast_server_output,n=5,species=None):
    """
    A slice of genuine NCBI blast output from the committed test data.

    Using real rows matters here: merge_and_annotate re-parses each hit's
    `title` with _parse_ncbi_line, so a hand-built title in the wrong format is
    silently dropped and the test ends up asserting nothing.
    """

    hits = ncbi_blast_server_output[0].iloc[:n].copy().reset_index(drop=True)
    hits["query"] = [f"count{i}" for i in range(len(hits))]

    if species is not None:
        hits["title"] = [f"ref|XP_00000{i}.1| some protein [{species}]"
                         for i in range(len(hits))]
        hits["hit_def"] = hits["title"]

    return hits


def test_df_from_seed_requires_a_blast_source(offline_seed):
    """
    With no NCBI db, no local db and no xml, there is nothing to blast against.
    """

    with pytest.raises(ValueError):
        df_from_seed(offline_seed,
                     ncbi_blast_db=None,
                     local_blast_db=None,
                     blast_xml=None)


def test_df_from_seed(mocker,offline_seed,ncbi_blast_server_output):
    """
    df_from_seed orchestrates blast -> merge -> ott -> nicknames. Check that the
    wiring between those steps does what it claims.
    """

    ncbi_blast = mocker.patch(
        "topiary.ncbi.ncbi_blast",
        return_value=[_real_blast_hits(ncbi_blast_server_output)])

    df, key_species, paralog_patterns, species_aware = df_from_seed(
        offline_seed,
        ncbi_blast_db="nr",
        hitlist_size=10)

    # BLAST was actually invoked against the database we asked for
    ncbi_blast.assert_called_once()
    assert ncbi_blast.call_args.kwargs["db"] == "nr"

    # Output is a valid topiary dataframe
    df = check.check_topiary_dataframe(df)
    for column in ["name","species","sequence","keep","uid","ott"]:
        assert column in df.columns

    # It holds the seed sequences *plus* the blast hits
    assert len(df) > len(key_species)
    assert "Homo sapiens" in list(df["species"])

    # uids are unique -- everything downstream keys on this
    assert len(set(df["uid"])) == len(df)

    # Seed sequences keep always_keep; new blast hits do not
    assert not df["always_keep"].all()
    assert df["always_keep"].any()

    # And the other three return values describe the seed
    assert "Homo sapiens" in key_species
    assert "LY96" in paralog_patterns
    assert isinstance(species_aware,bool)


def test_df_from_seed_drops_synthetic_sequences(mocker,offline_seed,ncbi_blast_server_output):
    """
    Synthetic constructs come back from NCBI with a real OTT that places them as
    an outgroup to all life, so they have to be dropped explicitly.
    """

    hits = _real_blast_hits(ncbi_blast_server_output,
                            n=2,
                            species="synthetic construct")

    mocker.patch("topiary.ncbi.ncbi_blast",return_value=[hits])

    df, _, _, _ = df_from_seed(offline_seed,ncbi_blast_db="nr")

    synth = df.loc[df["species"] == "synthetic construct",:]
    assert len(synth) > 0
    assert not synth["keep"].any()
