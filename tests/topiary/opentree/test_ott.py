
import pytest
import topiary
from topiary.opentree.ott import get_df_ott

import pandas as pd
import numpy as np

import re

@pytest.mark.network
def test_get_df_ott(test_dataframes):


    df = test_dataframes["good-df"]
    out_df = get_df_ott(df)
    assert out_df is not df

    tmp_df = df.drop(columns="ott")
    out_df = get_df_ott(tmp_df)
    assert np.array_equal(out_df.loc[:,"ott"],df.loc[:,"ott"])
    assert np.array_equal(out_df.loc[:,"species"],df.loc[:,"orig_species"])
    assert np.array_equal(np.ones(len(tmp_df),dtype=bool),out_df.keep)

    # make sure that it handles bad species names gracefully, setting keep to
    # False
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[tmp_df.index[0],"species"] = "Not a species"
    out_df = get_df_ott(tmp_df)
    assert pd.isnull(out_df.loc[out_df.index[0],"ott"])
    assert out_df.loc[out_df.index[0],"species"] == "Not a species"
    expected_keep = np.ones(len(out_df),dtype=bool)
    expected_keep[0] = False
    assert np.array_equal(expected_keep,out_df.keep)

    # make sure that it handles all bad species names gracefully, setting keep to
    # False
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[:,"species"] = "Not a species"
    out_df = get_df_ott(tmp_df)
    assert np.sum(pd.isnull(out_df.loc[:,"ott"])) == len(tmp_df)
    assert np.sum(out_df.loc[:,"species"] == "Not a species") == len(tmp_df)
    expected_keep = np.zeros(len(out_df),dtype=bool)
    assert np.array_equal(expected_keep,out_df.keep)

    tmp_df = df.drop(columns=["species","ott"])
    with pytest.raises(ValueError):
        out_df = get_df_ott(tmp_df)

    # make sure that it handles all bad species names gracefully, but keeps the
    # bad ott if we request it. 
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[:,"species"] = "Not a species"
    out_df = get_df_ott(tmp_df,keep_anyway=True)
    assert np.sum(pd.isnull(out_df.loc[:,"ott"])) == len(tmp_df)
    assert np.sum(out_df.loc[:,"species"] == "Not a species") == len(tmp_df)
    expected_keep = np.ones(len(out_df),dtype=bool)
    assert np.array_equal(expected_keep,out_df.keep)




# ---------------------------------------------------------------------------
# Hermetic tests. The network test above covers the happy path against the live
# Open Tree of Life database; these mock the two OTL calls so the branch logic
# around them -- keep flags, always_keep overrides, warnings -- is exercised in
# the default run rather than only when --run-network is given.
# ---------------------------------------------------------------------------

def _small_df():

    return pd.DataFrame({
        "name":["a","b","c"],
        "species":["Homo sapiens","Mus musculus","Not a species"],
        "sequence":["MAST","MASS","MASQ"],
        "keep":[True,True,True],
        "uid":["aaaaaaaaaa","bbbbbbbbbb","cccccccccc"],
    })


def _patch_otl(mocker,ott_list,resolvable,species_list=None):
    """
    Stand in for the two Open Tree of Life calls get_df_ott makes.
    """

    if species_list is None:
        species_list = ["Homo sapiens","Mus musculus","Not a species"]

    mocker.patch("topiary.opentree.ott.species_to_ott",
                 return_value=(ott_list,species_list,{}))
    mocker.patch("topiary.opentree.ott.ott_to_resolvable",
                 return_value=resolvable)


def test_get_df_ott_not_species_aware_skips_lookup(mocker):

    species_to_ott = mocker.patch("topiary.opentree.ott.species_to_ott")

    out = get_df_ott(_small_df(),species_aware=False)

    # No OTL traffic at all, and the columns are stubbed out
    species_to_ott.assert_not_called()
    assert out["ott"].isnull().all()
    assert out["resolvable"].all()

    # Nothing dropped
    assert out["keep"].all()

    # orig_species is still recorded
    assert list(out["orig_species"]) == list(_small_df()["species"])


def test_get_df_ott_validates_verbose():

    for bad in ["not_a_bool",1.5,None]:
        with pytest.raises(ValueError):
            get_df_ott(_small_df(),verbose=bad,species_aware=False)


def test_get_df_ott_drops_unrecognized_species(mocker,capsys):

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    out = get_df_ott(_small_df(),verbose=True)

    # The species with no ott is dropped, the others survive
    assert list(out["keep"]) == [True,True,False]
    assert pd.isnull(out.loc[out.index[2],"ott"])
    assert out.loc[out.index[0],"ott"] == "ott770315"

    # And the user is told why
    captured = capsys.readouterr()
    assert "Could not find OTT" in captured.out
    assert "Not a species" in captured.out


def test_get_df_ott_verbose_false_stays_quiet(mocker,capsys):

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    get_df_ott(_small_df(),verbose=False)

    captured = capsys.readouterr()
    assert "Could not find OTT" not in captured.out


def test_get_df_ott_drops_unresolvable_species(mocker,capsys):
    """
    A species can have a valid ott but still not be placeable on the synthetic
    tree (hybrids, for example). Those get dropped too, with a different
    message.
    """

    _patch_otl(mocker,
               ott_list=[770315,542509,12345],
               resolvable=[True,True,False])

    out = get_df_ott(_small_df(),verbose=True)

    assert list(out["keep"]) == [True,True,False]
    assert list(out["resolvable"]) == [True,True,False]

    captured = capsys.readouterr()
    assert "cannot be placed on tree" in captured.out


def test_get_df_ott_keep_anyway_leaves_keep_alone(mocker):

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    out = get_df_ott(_small_df(),keep_anyway=True,verbose=False)

    # Columns are still populated honestly...
    assert pd.isnull(out.loc[out.index[2],"ott"])

    # ...but nothing is dropped
    assert out["keep"].all()


def test_get_df_ott_overrides_always_keep_for_bad_species(mocker,capsys):
    """
    always_keep normally protects a sequence from quality filters, but a species
    that cannot be placed on the tree has to go regardless -- and the user is
    warned that always_keep was overridden.
    """

    df = _small_df()
    df["always_keep"] = [False,False,True]

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    out = get_df_ott(df,verbose=False)

    assert not bool(out.loc[out.index[2],"always_keep"])
    assert not bool(out.loc[out.index[2],"keep"])

    captured = capsys.readouterr()
    assert "always_keep" in captured.out


def test_get_df_ott_always_keep_respected_when_keep_anyway(mocker,capsys):

    df = _small_df()
    df["always_keep"] = [False,False,True]

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    out = get_df_ott(df,keep_anyway=True,verbose=False)

    # keep_anyway suppresses both the drop and the warning
    assert out["keep"].all()
    captured = capsys.readouterr()
    assert "always_keep" not in captured.out


def test_get_df_ott_returns_a_copy(mocker):

    _patch_otl(mocker,ott_list=[770315,542509,None],resolvable=[True,True])

    df = _small_df()
    out = get_df_ott(df,verbose=False)

    assert out is not df

    # The input is untouched
    assert "ott" not in df.columns
    assert df["keep"].all()
