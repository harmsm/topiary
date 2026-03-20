
import pytest

import topiary
from topiary.quality import taxonomic as tx

import ete4 as ete
import numpy as np
import pandas as pd

import random, os

def test__prep_species_tree(test_dataframes):

    df = test_dataframes["good-df"].copy()
    df["recip_paralog"] = "LY96"

    df, T = tx._prep_species_tree(df,paralog_column="recip_paralog")
    assert type(T) is ete.Tree
    for leaf in T.leaves():
        assert np.array_equal(list(leaf.get_prop("paralogs").keys()),["LY96"])
        assert len(leaf.get_prop("paralogs")["LY96"]) == 1
        assert leaf.get_prop("uid")[0] == leaf.get_prop("paralogs")["LY96"][0] # Make sure we're not mixing up uid

    with pytest.raises(ValueError):
        df, T = tx._prep_species_tree(df,paralog_column="not_a_column")

    # Test ott generation if missing
    df_no_ott = test_dataframes["good-df"].copy()
    df_no_ott = df_no_ott.drop(columns=["ott"])
    df_no_ott["recip_paralog"] = "LY96"
    df_out, T_out = tx._prep_species_tree(df_no_ott, paralog_column="recip_paralog")
    assert "ott" in df_out.columns

    df = test_dataframes["good-df"].copy()
    second_df = test_dataframes["good-df"].copy()

    df = pd.concat((df,second_df),ignore_index=True)
    df.uid = topiary._private.generate_uid(len(df))
    df["recip_paralog"] = "LY96"
    df.loc[0:4,"recip_paralog"] = "LY86"

    df, T = tx._prep_species_tree(df,paralog_column="recip_paralog")
    for leaf in T.leaves():
        paralog_keys = list(leaf.get_prop("paralogs").keys())
        paralog_keys.sort()
        assert np.array_equal(paralog_keys,["LY86","LY96"])
        assert len(leaf.get_prop("paralogs")["LY96"]) == 1
        assert len(leaf.get_prop("paralogs")["LY86"]) == 1

def test__even_paralog_budgeting():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make tree with X and Y paralogs at the tip, each of which are seen once
    T = ete.Tree(some_tree)
    for leaf in T.leaves():
        leaf.add_prop("paralogs", {})
        leaf.get_prop("paralogs")["X"] = [topiary._private.generate_uid(1)]
        leaf.get_prop("paralogs")["Y"] = [topiary._private.generate_uid(1)]
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=18)
    assert budget["X"] == 7
    assert budget["Y"] == 7

def test__weighted_paralog_budgeting():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    T = ete.Tree(some_tree)
    for leaf in T.leaves():
        leaf.add_prop("paralogs", {})
        leaf.get_prop("paralogs")["X"] = [topiary._private.generate_uid(1)]
        leaf.get_prop("paralogs")["Y"] = [topiary._private.generate_uid(1)]
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=18)
    assert budget["X"] == 7
    assert budget["Y"] == 7

def test__finalize_paralog_budget():

    paralog_budget = {'Y': 10, 'X': 11}
    paralog_counts = {'Y': 14, 'X': 7}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 7
    assert final_budget["Y"] == 14

def test__get_sequence_budgets():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    T = ete.Tree(some_tree)
    for leaf in T.leaves():
        leaf.add_prop("sequences",("X",))
    T_bak = T.copy()

    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,7)
    assert T.get_prop("budget") == 7

def test__even_merge_blocks():
    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    T = ete.Tree(some_tree)
    for leaf in T.leaves():
        leaf.add_prop("sequences",(topiary._private.uid.generate_uid(1),))
        leaf.add_prop("num_seq", 1)
    
    # Manually set num_seq and budget for nodes
    for node in T.traverse():
        if not node.is_leaf:
            node.add_prop("num_seq", len(list(node.leaves())))
    
    # Merge block size 3
    blocks = tx._even_merge_blocks(T, 3)
    assert len(blocks) > 0

def test_get_merge_blocks(test_dataframes):
    df = test_dataframes["good-df"].copy()
    df["recip_paralog"] = "P1"
    
    # Test dummy_merge_blocks
    blocks = tx.get_merge_blocks(df, target_seq_number=3, dummy_merge_blocks=True)
    assert None in blocks
    
    # Test ValueError for invalid column
    with pytest.raises(ValueError):
        tx.get_merge_blocks(df, target_seq_number=3, paralog_column="missing")
    
    # Test standard run
    blocks = tx.get_merge_blocks(df, target_seq_number=2)
    assert "P1" in blocks

    # Test weighted_paralog_split
    blocks = tx.get_merge_blocks(df, target_seq_number=2, weighted_paralog_split=True)
    assert "P1" in blocks

    # Test target_merge_block_size
    blocks = tx.get_merge_blocks(df, target_seq_number=2, target_merge_block_size=1)
    assert "P1" in blocks

    # Test RuntimeError in _taxonomic_merge_blocks
    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    T = ete.Tree(some_tree)
    for leaf in T.leaves():
        leaf.add_prop("sequences",(topiary._private.uid.generate_uid(1),))
    
    # Test budgeting for 1 sequence
    T = tx._get_sequence_budgets(T, 1)
