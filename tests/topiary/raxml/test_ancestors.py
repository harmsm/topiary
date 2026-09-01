import pytest

import topiary
from topiary.raxml.ancestors import _make_ancestor_summary_trees
from topiary.raxml.ancestors import _parse_raxml_anc_output
from topiary.raxml.ancestors import _get_bad_columns
from topiary.raxml.ancestors import generate_ancestors
from topiary.raxml import RAXML_BINARY
from topiary._private.interface import WrappedFunctionException

import numpy as np

import os
import copy
import json

def test__make_ancestor_summary_trees(tmpdir):
    """
    Writes two newick files from a labeled raxml tree: one with internal nodes
    renamed to ancestor names, one with them replaced by average posterior
    probability.
    """

    from ete4 import Tree

    os.chdir(tmpdir)

    # Internal nodes named the way raxml labels them
    labeled = "((A:0.1,B:0.1)Node1:0.5,(C:0.1,D:0.1)Node2:0.5)Node3:1.0;"
    with open("labeled.newick","w") as f:
        f.write(labeled)

    avg_pp_dict = {"anc1":0.95,"anc2":0.85,"anc3":0.75}

    _make_ancestor_summary_trees(None,avg_pp_dict,"labeled.newick")

    assert os.path.isfile("tree_anc-label.newick")
    assert os.path.isfile("tree_anc-pp.newick")

    # Label tree: NodeN becomes ancN
    t_label = Tree(open("tree_anc-label.newick").read(),parser=1)
    internal_names = set([n.name for n in t_label.traverse() if not n.is_leaf])
    assert "anc1" in internal_names
    assert "anc2" in internal_names
    assert "Node1" not in internal_names

    # Leaves are untouched
    assert set([n.name for n in t_label.leaves()]) == {"A","B","C","D"}

    # pp tree: internal nodes carry the posterior probabilities
    t_pp = Tree(open("tree_anc-pp.newick").read(),parser=1)
    pp_names = set([n.name for n in t_pp.traverse() if not n.is_leaf])
    assert "0.95" in pp_names
    assert "0.85" in pp_names


def test__parse_raxml_anc_output(tmpdir,small_phylo):
    """
    Parse a real raxml marginal-ancestor run out of the committed test data and
    check the human-readable outputs it is supposed to produce.
    """

    os.chdir(tmpdir)

    df = topiary.read_dataframe(small_phylo["02_gene-tree-ancestors/input/dataframe.csv"])
    anc_prob_file = small_phylo["02_gene-tree-ancestors/working/00_inference/alignment.phy.raxml.ancestralProbs"]
    alignment_file = small_phylo["02_gene-tree-ancestors/working/00_inference/alignment.phy"]
    tree_with_labels = small_phylo["02_gene-tree-ancestors/working/00_inference/alignment.phy.raxml.ancestralTree"]

    # The caller (generate_ancestors, via Supervisor) creates the output
    # directory before handing it to the parser.
    os.mkdir("ancestors")

    _parse_raxml_anc_output(df,
                            anc_prob_file,
                            alignment_file,
                            tree_with_labels,
                            run_directory="ancestors")

    # The parser writes a directory of human-readable ancestor output
    assert os.path.isdir("ancestors")

    written = os.listdir("ancestors")

    # A fasta of ancestral sequences and a csv of per-site probabilities
    assert any(f.endswith(".fasta") for f in written)
    assert any(f.endswith(".csv") for f in written)

    # Plus the two summary trees
    assert "tree_anc-label.newick" in written
    assert "tree_anc-pp.newick" in written

    # The fasta must actually contain ancestor sequences, not be empty
    fasta = [f for f in written if f.endswith(".fasta")][0]
    with open(os.path.join("ancestors",fasta)) as f:
        contents = f.read()
    assert contents.startswith(">")
    assert len(contents.strip().split("\n")) >= 2

def test__get_bad_columns(tmpdir):

    os.chdir(tmpdir)

    bad_file = ["",""]
    bad_file.extend(["seq1","AXAAAXXXA-AAA---X-XXX---"])
    bad_file.extend(["seq2","AAXAXAXXAA-A-A--XX-X-X--"])
    bad_file.extend(["seq3","AAAXXXAXAAA---A-XXX---X-"])
    expected_bad = np.array([int(x) for x in list("000000010000000111111111")],dtype=bool)
    expected_bad = np.arange(len(expected_bad),dtype=int)[expected_bad]

    f = open("bad-file.phy","w")
    f.write("\n".join(bad_file))
    f.close()

    observed_bad = _get_bad_columns("bad-file.phy")
    assert np.array_equal(observed_bad,expected_bad)


@pytest.mark.run_raxml
def test_generate_ancestors(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]
    gene_tree = tiny_phylo["final-output/gene-tree.newick"]
    species_tree = tiny_phylo["initial-input/species-tree.newick"]
    reconciled_tree = tiny_phylo["final-output/reconciled-tree.newick"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"prev_calculation":None,
                       "df":df,
                       "model":"JTT",
                       "gene_tree":gene_tree,
                       "alt_cutoff":0.25,
                       "calc_dir":"ancestors",
                       "overwrite":False,
                       "num_threads":1,
                       "raxml_binary":RAXML_BINARY}

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["calc_dir"] = "test0"

    generate_ancestors(**kwargs)

    output = os.path.join("test0","output")

    expected = ["gene-tree_anc-label.newick",
                "gene-tree_anc-pp.newick",
                "summary-tree.pdf",
                "dataframe.csv"]

    expected.append(os.path.join("gene-tree_ancestors","ancestors.fasta"))
    expected.append(os.path.join("gene-tree_ancestors","ancestor-data.csv"))
    for i in range(6):
        expected.append(os.path.join("gene-tree_ancestors",f"anc{i+1}.pdf"))

    for e in expected:
        assert os.path.isfile(os.path.join(output,e))

    f = open(os.path.join("test0","run_parameters.json"))
    run_params = json.load(f)
    f.close()
    assert run_params["model"] == "JTT"
    assert run_params["alt_cutoff"] == 0.25


    # --------------------------------------------------------------------------
    # Make sure reconciled tree takes precendence over gene tree

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["reconciled_tree"] = reconciled_tree
    kwargs["model"] = "LG"
    kwargs["alt_cutoff"] = 0.20
    kwargs["calc_dir"] = "test1"

    generate_ancestors(**kwargs)

    output = os.path.join("test1","output")

    expected = ["reconciled-tree_anc-label.newick",
                "reconciled-tree_anc-pp.newick",
                "summary-tree.pdf",
                "dataframe.csv"]

    expected.append(os.path.join("reconciled-tree_ancestors","ancestors.fasta"))
    expected.append(os.path.join("reconciled-tree_ancestors","ancestor-data.csv"))
    for i in range(6):
        expected.append(os.path.join("reconciled-tree_ancestors",f"anc{i+1}.pdf"))

    for e in expected:
        assert os.path.isfile(os.path.join(output,e))

    f = open(os.path.join("test1","run_parameters.json"))
    run_params = json.load(f)
    f.close()
    assert run_params["model"] == "LG"
    assert run_params["alt_cutoff"] == 0.20

    os.chdir("..")

    # --------------------------------------------------------------------------
    # make sure it dies when neither gene tree nor reconciled passed in

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["calc_dir"] = "test2"
    kwargs["gene_tree"] = None
    kwargs["reconciled_tree"] = None

    with pytest.raises(ValueError):
        generate_ancestors(**kwargs)
