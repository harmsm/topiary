import pytest
import topiary
from topiary.opentree.util import _validate_ott_or_species
from topiary.opentree.util import ott_to_mrca
from topiary.opentree.util import ott_to_resolvable
from topiary.opentree.util import species_to_ott
from topiary.opentree.util import _tnrs_query
from topiary.opentree.util import _parse_tnrs_matches
from topiary.opentree.util import ott_to_species_tree
from topiary.opentree.util import tree_to_taxa_order
from topiary.opentree.util import sort_df_by_taxa

import ete4 as ete

import pandas as pd
import numpy as np

import re
import string
import warnings

@pytest.mark.network
def test__validate_ott_or_species():

    # Nothing in
    with pytest.raises(ValueError):
        _validate_ott_or_species(None,None)

    # Both in
    with pytest.raises(ValueError):
        _validate_ott_or_species([],[])

    # Bad ott
    with pytest.raises(ValueError):
        _validate_ott_or_species(ott_list=["Homo sapiens","Gallus gallus"])

    # Bad species
    with pytest.raises(ValueError):
        _validate_ott_or_species(species_list=[111,222])

    # Bad species
    with pytest.raises(ValueError):
        _validate_ott_or_species(species_list=1.0)

    # Bad ott
    with pytest.raises(ValueError):
        _validate_ott_or_species(ott_list=["ott"])

    # Bad ott
    with pytest.raises(ValueError):
        _validate_ott_or_species(ott_list=1.0)

    # Good ott
    ott_list = _validate_ott_or_species(ott_list=[111,222])
    assert len(ott_list) == 2
    assert ott_list[0] == 111
    assert ott_list[1] == 222

    # Good ott
    ott_list = _validate_ott_or_species(ott_list=["ott111",222])
    assert len(ott_list) == 2
    assert ott_list[0] == 111
    assert ott_list[1] == 222

    # Good species
    ott_list = _validate_ott_or_species(species_list=["Homo sapiens"])
    assert len(ott_list) == 1
    assert ott_list[0] == 770315

    # empty ott
    ott_list = _validate_ott_or_species(ott_list=[])
    assert len(ott_list) == 0

    # empty species
    ott_list = _validate_ott_or_species(species_list=[])
    assert len(ott_list) == 0


@pytest.mark.network
def test_species_to_ott():

    input_species_list = ["Homo sapiens",
                          "Thermus thermophilus",
                          "Methanococcus voltae",
                          "Saccharomyces cerevisiae"]

    expected_ott = [770315,276534,565131,356221]
    expected_taxid = [9606,274,2188,4932]

    ott_list, species_list, results = species_to_ott(input_species_list)
    assert len(ott_list) == 4
    assert len(species_list) == 4
    assert len(results) == 4
    assert np.array_equal(ott_list,expected_ott)
    assert np.array_equal(species_list,input_species_list)
    assert issubclass(type(results),dict)
    for i, s in enumerate(input_species_list):
        assert results[s]["num_matches"] == 1
        assert results[s]["taxid"] == expected_taxid[i]

    # Make sure it checks for inputs correctly
    bad_inputs = [None,[1.5,2.5],1,{1:"test"},list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            species_to_ott(b)

    # Make sure it handles an empty input list gracefully
    ott_list, species_list, results = species_to_ott([])
    assert len(ott_list) == 0
    assert len(species_list) == 0

    # An unresolvable name warns the user and yields ott None. The warning is
    # part of the contract, so assert on it rather than letting it escape.
    with pytest.warns(UserWarning,match="could not be assigned"):
        ott_list, species_list, results = species_to_ott(["Not really a species"])
    assert len(ott_list) == 1
    assert len(species_list) == 1
    assert ott_list[0] is None
    assert species_list[0] == "Not really a species"\

    # Two ambiguous cases. Gekko has a synonym an should be resolved. Capra
    # has a (wacky) bacterial species name and should not resolve. If this test
    # starts failing, opentreeoflife may finally have fixed problem where Capra
    # resolves to a bacterium...
    ott_list, species_list, results = species_to_ott(["Gekko japonicus",
                                                      "Capra hircus",
                                                      "Escherichia coli"])
    assert len(ott_list) == 3
    assert ott_list[0] == 212506
    assert ott_list[1] == 19017
    assert ott_list[2] == 474506

    # This should work because we're inside cow-like things. OTT should realize
    # we mean a goat, not a bacterium.
    ott_list, species_list, results = species_to_ott(["Bos taurus",
                                                      "Capra hircus"])
    assert len(ott_list) == 2
    assert ott_list[0] == 490099
    assert ott_list[1] == 19017

    # This is a fuzzy match and should fail -- warning included
    with pytest.warns(UserWarning,match="could not be assigned"):
        ott_list, species_list, results = species_to_ott(["Neosciurus carolinensis"])
    assert results["Neosciurus carolinensis"]["msg"].startswith("No exact match")


def _fake_tnrs_response(name_to_ott):
    """
    Build a minimal stand-in for the object OT.tnrs_match returns. Takes a
    dictionary of ott ids keyed to query name; an ott of None means the query
    matched nothing.
    """

    class FakeResponse:
        pass

    results = []
    for name, ott_id in name_to_ott.items():
        if ott_id is None:
            matches = []
        else:
            matches = [{"is_approximate_match":False,
                        "is_synonym":False,
                        "matched_name":name,
                        "taxon":{"ott_id":ott_id,
                                 "name":name,
                                 "tax_sources":[]}}]
        results.append({"name":name,"matches":matches})

    out = FakeResponse()
    out.response_dict = {"results":results}

    return out


def test__tnrs_query(monkeypatch):

    # An empty query should not hit the server at all
    def explode(*args,**kwargs):
        raise AssertionError("should not have queried opentree")
    monkeypatch.setattr(topiary.opentree.util.OT,"tnrs_match",explode)
    assert _tnrs_query([]) == {}

    # Results should be keyed by the name each result carries, not by the
    # position it came back in. Hand back results in an order that has nothing
    # to do with the query order.
    sent = {}
    def fake_tnrs_match(names,context_name=None,do_approximate_matching=False,
                        include_suppressed=False):
        sent["names"] = names
        sent["context_name"] = context_name
        return _fake_tnrs_response({"Homo sapiens":770315,
                                    "Octopus sinensis":110468,
                                    "Hypanus sabinus":351766})

    monkeypatch.setattr(topiary.opentree.util.OT,"tnrs_match",fake_tnrs_match)

    query = ["Octopus sinensis","Hypanus sabinus","Homo sapiens"]
    matches = _tnrs_query(query)

    # We should have asked for the widest possible context
    assert sent["names"] == query
    assert sent["context_name"] == topiary.opentree.util._TNRS_CONTEXT

    for s, expected_ott in [("Octopus sinensis",110468),
                            ("Hypanus sabinus",351766),
                            ("Homo sapiens",770315)]:
        assert matches[s][0]["taxon"]["ott_id"] == expected_ott


def test__parse_tnrs_matches():

    # No match
    result = _parse_tnrs_matches("Not really a species",[])
    assert result["matched"] == False
    assert result["num_matches"] == 0
    assert result["ott_id"] is None
    assert result["ott_name"] is None
    assert result["taxid"] is None
    assert result["msg"].startswith("No match")

    # Single, exact match
    matches = _fake_tnrs_response({"Homo sapiens":770315}).response_dict["results"][0]["matches"]
    matches[0]["taxon"]["tax_sources"] = ["ncbi:9606","gbif:2436436"]
    result = _parse_tnrs_matches("Homo sapiens",matches)
    assert result["matched"] == True
    assert result["num_matches"] == 1
    assert result["msg"] == "success"
    assert result["ott_id"] == 770315
    assert result["ott_name"] == "Homo sapiens"
    assert result["taxid"] == 9606

    # Multiple, non-synonymous matches are ambiguous
    matches = _fake_tnrs_response({"Drosophila":1}).response_dict["results"][0]["matches"]
    matches[0]["taxon"]["tax_sources"] = ["ncbi:7215"]
    matches.append({"is_approximate_match":False,
                    "is_synonym":False,
                    "matched_name":"Drosophila",
                    "taxon":{"ott_id":2,
                             "name":"Drosophila",
                             "tax_sources":["ncbi:2081351"]}})
    result = _parse_tnrs_matches("Drosophila",matches)
    assert result["num_matches"] == 2
    assert result["msg"] == "multiple matches"
    assert result["ott_id"] == 1

    # An approximate match is not a match
    matches = _fake_tnrs_response({"Neosciurus carolinensis":1}).response_dict["results"][0]["matches"]
    matches[0]["is_approximate_match"] = True
    matches[0]["matched_name"] = "Sciurus carolinensis"
    matches[0]["taxon"]["tax_sources"] = ["ncbi:1"]
    matches.append({"is_approximate_match":True,
                    "is_synonym":False,
                    "matched_name":"Nesticus carolinensis",
                    "taxon":{"ott_id":2,
                             "name":"Nesticus carolinensis",
                             "tax_sources":["ncbi:2"]}})
    result = _parse_tnrs_matches("Neosciurus carolinensis",matches)
    assert result["matched"] == False
    assert result["msg"].startswith("No exact match")


@pytest.mark.network
def test_species_to_ott_batch_context():
    """
    opentree's TNRS infers a taxonomic context from the batch of names sent to
    it. Historically, mixing distant clades in one query narrowed that context
    and made the out-of-context species come back unmatched, even though those
    same species matched fine when queried alone or in a smaller batch.
    """

    expected = {"Hypanus sabinus":351766,
                "Octopus sinensis":110468,
                "Homo sapiens":770315}

    queries = [["Hypanus sabinus"],
               ["Octopus sinensis"],
               ["Hypanus sabinus","Octopus sinensis"],
               ["Octopus sinensis","Homo sapiens","Hypanus sabinus"]]

    for query in queries:
        ott_list, species_list, results = species_to_ott(query)
        assert len(ott_list) == len(query)
        for i, s in enumerate(query):
            assert ott_list[i] == expected[s]
            assert results[s]["num_matches"] == 1


@pytest.mark.network
def test_species_to_ott_warns_on_unmatched():

    # A species that cannot be matched should not vanish silently
    with pytest.warns(UserWarning,match="Not really a species"):
        ott_list, species_list, results = species_to_ott(["Homo sapiens",
                                                          "Not really a species"])
    assert ott_list[0] == 770315
    assert ott_list[1] is None

    # ... and one that can be matched should not warn
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ott_list, species_list, results = species_to_ott(["Homo sapiens"])
    assert ott_list[0] == 770315


@pytest.mark.network
def test_species_to_ott_strains():
    # Test strain stripping
    input_species_list = ["Escherichia coli (strain K12)",
                          "Escherichia coli strain K12",
                          "Escherichia coli str. K12",
                          "Vibrio cholerae O1 biovar El Tor str. N16961"]
    
    # All E. coli should map to Escherichia coli (ott 474506)
    # The Vibrio should map to Vibrio cholerae (ott 671232)
    
    ott_list, species_list, results = species_to_ott(input_species_list)
    
    # Check E. coli results
    for i in range(3):
        assert ott_list[i] == 474506
        assert species_list[i] == "Escherichia coli"
        assert results[input_species_list[i]]["matched"] == True
        
    # Check Vibrio result
    assert ott_list[3] == 77943
    # Note: this might be 'Vibrio cholerae O1 biovar El Tor' since only the
    # 'str. N16961' was stripped.
    assert "Vibrio cholerae" in species_list[3]


@pytest.mark.network
def test_ott_to_species_tree():

    input_ott = [770315,276534,565131,356221]

    T, results = ott_to_species_tree(input_ott)
    assert issubclass(type(T),ete.Tree)
    assert issubclass(type(results),dict)
    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()

    assert np.array_equal(seen,sent_in)
    assert len(results["not_resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0

    off_tree = []
    for n in T.leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in)

     # Send in a bad ott
    input_ott = [770315,276534,565131,356221,9999999999999999999999]
    T, results = ott_to_species_tree(input_ott)

    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()
    assert np.array_equal(seen,sent_in[:-1])

    off_tree = []
    for n in T.leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in[:-1])

    assert len(results["not_resolved"]) == 1
    assert results["not_resolved"][0] == 9999999999999999999999
    assert len(results["unknown_ids"]) == 1
    assert results["unknown_ids"][0] == 9999999999999999999999
    assert len(results["not_monophyletic"]) == 0

    # Make sure it handles empty list gracefully
    T, results = ott_to_species_tree([])
    assert T is None
    assert len(results["resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0
    assert len(results["not_resolved"]) == 0

     # Send in only one ott
    input_ott = [770315]
    T, results = ott_to_species_tree(input_ott)

    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()
    assert np.array_equal(seen,sent_in)

    off_tree = []
    for n in T.leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in)

    assert len(results["not_resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0

     # Send in only one ott, but it's bad
    input_ott = [9999999999999999999999]
    T, results = ott_to_species_tree(input_ott)
    assert T is None
    assert len(results["resolved"]) == 0
    assert len(results["not_resolved"]) == 1
    assert results["not_resolved"][0] == 9999999999999999999999
    assert len(results["unknown_ids"]) == 1
    assert results["unknown_ids"][0] == 9999999999999999999999
    assert len(results["not_monophyletic"]) == 0

@pytest.mark.network
def test_ott_to_resolvable():

    some_good_ott = [770315,276534,565131,356221]

    # hybrid that is not on synthetic tree
    bad_ott = [4942641]

    resolvable = ott_to_resolvable(some_good_ott)
    assert np.array_equal(np.ones(len(some_good_ott),dtype=bool),resolvable)

    resolvable = ott_to_resolvable(bad_ott)
    assert np.array_equal(np.zeros(len(bad_ott),dtype=bool),resolvable)

    both_together = some_good_ott[:]
    both_together.extend(bad_ott)
    resolvable = ott_to_resolvable(both_together)
    expected = np.array([1,1,1,1,0],dtype=bool)
    assert np.array_equal(resolvable,expected)

    replicated = [770315,770315,770315,770315,770315,4942641,4942641,4942641,4942641]
    expected = np.array([1,1,1,1,1,0,0,0,0],dtype=bool)
    resolvable = ott_to_resolvable(replicated)
    assert np.array_equal(expected,resolvable)

    # Make sure it handes one good taxa well
    resolvable = ott_to_resolvable([770315])
    assert len(resolvable) == 1
    assert resolvable[0] == True

    # Make sure it handes one bad taxa well
    resolvable = ott_to_resolvable([4942641])
    assert len(resolvable) == 1
    assert resolvable[0] == False

    # Make sure it handles single empty list
    resolvable = ott_to_resolvable([])
    assert len(resolvable) == 0

    # Make sure it checks for inputs correctly
    bad_inputs = [None,1,list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            ott_to_resolvable(b)


@pytest.mark.network
def test_ott_to_mrca():

    vert = [770315, 153563, 1005914]           #["Homo sapiens","Gallus gallus","Danio rerio"]
    amniotes = [770315, 153563]                #["Homo sapiens","Gallus gallus"]
    tetrapods = [770315, 153563, 465096]       #["Homo sapiens","Gallus gallus","Xenopus laevis"]
    placental_mammals = [770315, 542509]       #["Homo sapiens","Mus musculus"]
    therian_mammals = [770315, 542509, 122362] #["Homo sapiens","Mus musculus","Monodelphis domestica"]
    plants = [309263,1075713]                  #["Arabidopsis thaliana","Pseudotsuga menziesii"]
    eukaryotes = [770315, 356221, 1075713]     #["Homo sapiens","Saccharomyces cerevisiae","Pseudotsuga menziesii"]
    human_yeast = [770315, 356221]             #["Homo sapiens","Saccharomyces cerevisiae"]
    archaea = [7000614,996421]                 #["Mancarchaeum acidiphilum","Mendosicutes"]
    not_bacteria = [770315,7000614,996421]     #["Homo sapiens","Mancarchaeum acidiphilum","Mendosicutes"]
    all_life = [770315, 474506, 7000614]       #["Homo sapiens","Escherichia coli","Mancarchaeum acidiphilum"]
    bacteria = [1090496, 474506]               #["Staphylococcus aureus","Escherichia coli"]

    # Get with default args
    out = ott_to_mrca(vert)
    assert out["ott_name"] == 'Euteleostomi'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 117571
    assert out["is_microbial"] == False

    # Move up four ranks from Euteleostomi
    out = ott_to_mrca(vert,move_up_by=4)
    assert out["ott_name"] == 'Craniata'
    assert out["ott_rank"] == 'subphylum'
    assert out["taxid"] == 89593
    assert out["is_microbial"] == False

    # Move waaaaay up but do not allow all life (end up on Eukaryota)
    out = ott_to_mrca(vert,move_up_by=10000)
    assert out["ott_name"] == 'Eukaryota'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2759
    assert out["is_microbial"] == False

    # Move waaaaay up, allowing all life
    out = ott_to_mrca(vert,avoid_all_life=False,move_up_by=10000)
    assert out["ott_name"] == 'life'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 1
    assert out["is_microbial"] == False

    # Try to avoid all life, but we can't because inputs are from all three
    # domains.
    out = ott_to_mrca(all_life,avoid_all_life=True)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567
    assert out["is_microbial"] == False

    # Spot check various taxonomic groups
    out = ott_to_mrca(amniotes)
    assert out["ott_name"] == 'Amniota'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 32524
    assert out["is_microbial"] == False

    out = ott_to_mrca(tetrapods)
    assert out["ott_name"] == 'Tetrapoda'
    assert out["ott_rank"] == 'superclass'
    assert out["taxid"] == 32523
    assert out["is_microbial"] == False

    out = ott_to_mrca(placental_mammals)
    assert out["ott_name"] == 'Euarchontoglires'
    assert out["ott_rank"] == 'superorder'
    assert out["taxid"] == 314146
    assert out["is_microbial"] == False

    out = ott_to_mrca(therian_mammals)
    assert out["ott_name"] == 'Theria'
    assert out["ott_rank"] == 'subclass'
    assert out["taxid"] == 32525
    assert out["is_microbial"] == False

    out = ott_to_mrca(human_yeast)
    assert out["ott_name"] == 'Opisthokonta'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 33154
    assert out["is_microbial"] == False

    out = ott_to_mrca(plants)
    assert out["ott_name"] == 'Spermatophyta'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 58024
    assert out["is_microbial"] == False

    out = ott_to_mrca(eukaryotes)
    assert out["ott_name"] == 'Eukaryota'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2759
    assert out["is_microbial"] == False

    out = ott_to_mrca(archaea)
    assert out["ott_name"] == 'Archaea'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2157
    assert out["is_microbial"] == True

    out = ott_to_mrca(bacteria)
    assert out["ott_name"] == 'Bacteria'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2
    assert out["is_microbial"] == True

    out = ott_to_mrca(not_bacteria)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567
    assert out["is_microbial"] == False

    out = ott_to_mrca(all_life)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567
    assert out["is_microbial"] == False

    # Make sure it checks for inputs correctly
    bad_inputs = [None,1,list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            ott_to_mrca(b)

    # Pass in bad values for move_up_by and avoid_all_life
    with pytest.raises(ValueError):
        ott_to_mrca(vert,move_up_by=-1)

    with pytest.raises(ValueError):
        ott_to_mrca(vert,move_up_by="stupid")

    with pytest.raises(ValueError):
        ott_to_mrca(vert,avoid_all_life="True")

    # Send in empty list
    out = ott_to_mrca([])
    assert out["ott_name"] == 'life'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 1

    # Send in only bad ott
    with pytest.raises(ValueError):
        out = ott_to_mrca([99999999999999])

    # Send in placentl mammals with a bad ott
    to_test = placental_mammals[:]
    to_test.append(99999999999999)
    out = ott_to_mrca(to_test)
    assert out["ott_name"] == 'Euarchontoglires'
    assert out["ott_rank"] == 'superorder'
    assert out["taxid"] == 314146
    assert out["is_microbial"] == False

def test_tree_to_taxa_order():

    T = ete.Tree("((((A,B),(C,D)),Q),(H,(E,F)));")
    out_order = tree_to_taxa_order(T,ref_name="H")
    assert np.array_equal(out_order[:4],["H","E","F","Q"])

    out_order = tree_to_taxa_order(T,ref_name="A")
    assert np.array_equal(out_order,["A","B","C","D","Q","H","E","F"])

    out_order = tree_to_taxa_order(T,ref_name="B")
    assert np.array_equal(out_order,["B","A","C","D","Q","H","E","F"])

    out_order = tree_to_taxa_order(T,ref_name="F")
    assert np.array_equal(out_order[:4],["F","E","H","Q"])

    # ref_name not in tree
    out_order = tree_to_taxa_order(T,ref_name="X")
    assert isinstance(out_order,list)
    assert len(out_order) == 8

    # No ref_name given
    out_order = tree_to_taxa_order(T)
    assert isinstance(out_order,list)
    assert len(out_order) == 8

@pytest.mark.network
def test_sort_df_by_taxa(for_real_inference):

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])

    # Make sure it sorts on first key_species ott then recip_paralog
    new_df = sort_df_by_taxa(df)
    assert new_df.loc[new_df.index[0],"ott"] == "ott770315"
    expected = ["LY86","LY86","LY86","LY86","LY86","LY86","LY86","LY86",
                "LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96",
                "unassigned"]
    assert np.array_equal(new_df.loc[:,"recip_paralog"],expected)

    # Make sure it sorts on given ott then recip_paralog
    new_df = sort_df_by_taxa(df,ref_ott="ott490109")
    assert new_df.loc[new_df.index[0],"ott"] == "ott490109"
    expected = ["LY86","LY86","LY86","LY86","LY86","LY86","LY86","LY86",
                "LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96","LY96",
                "unassigned"]
    assert np.array_equal(new_df.loc[:,"recip_paralog"],expected)

    # Drop recip_paralog column; should now sort on nickname
    input_df = df.copy()
    input_df = input_df.drop(columns=["recip_paralog"])
    values = [str(x).zfill(3) for x in range(len(input_df))]
    np.random.shuffle(values)
    input_df.loc[:,"nickname"] = values
    new_df = sort_df_by_taxa(input_df,ref_ott="ott770315")
    assert np.array_equal(new_df.loc[:,"nickname"],[str(x).zfill(3) for x in range(len(input_df))])

    # Drop recip_paralog and nickname column; should now sort on name
    input_df = df.copy()
    input_df = input_df.drop(columns=["recip_paralog","nickname"])
    values = list(string.ascii_lowercase[:len(input_df)])
    np.random.shuffle(values)
    input_df.loc[:,"name"] = values
    new_df = sort_df_by_taxa(input_df,ref_ott="ott770315")
    assert np.array_equal(new_df.loc[:,"name"],list(string.ascii_lowercase[:len(input_df)]))

    # Short on custom paralog column
    input_df = df.copy()
    values = [str(x).zfill(3) for x in range(len(input_df))]
    np.random.shuffle(values)
    input_df["rocket"] = values
    new_df = sort_df_by_taxa(input_df,paralog_column="rocket")
    assert np.array_equal(new_df.loc[:,"rocket"],[str(x).zfill(3) for x in range(len(input_df))])

    # Make sure it's paying attention to only_keepers
    input_df = df.copy()
    input_df.loc[input_df.index[5:],"keep"] = False
    new_df = sort_df_by_taxa(input_df)
    assert len(new_df) == len(input_df)

    input_df = df.copy()
    input_df.loc[input_df.index[5:],"keep"] = False
    new_df = sort_df_by_taxa(input_df,only_keepers=True)
    assert len(new_df) == 5

    input_df = df.copy()
    input_df.loc[input_df.index[5:],"keep"] = False
    new_df = sort_df_by_taxa(input_df,only_keepers=False)
    assert len(new_df) == len(input_df)

    # Arg checking
    bad_paralog_columns = ["A",1.5,list]
    for b in bad_paralog_columns:
        print("passing",b)
        with pytest.raises(ValueError):
            new_df = sort_df_by_taxa(df,paralog_column=b)

    bad_ott = [770315,"Homo sapiens",1.1]
    for b in bad_ott:
        print("passing",b)
        with pytest.raises(ValueError):
            new_df = sort_df_by_taxa(df,ref_ott=b)

    bad_keep = ["Homo sapiens",1.1,None]
    for b in bad_keep:
        print("passing",b)
        with pytest.raises(ValueError):
            new_df = sort_df_by_taxa(df,only_keepers=b)
