import pytest
from topiary.reports.cards.duplications import create_duplications_card, _check_duplication
from topiary._private.supervisor import Supervisor
from topiary.io.tree import load_trees
import os
import pandas as pd

def test_duplications(tiny_phylo):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["03_ancestors"]
    sv = Supervisor(calc_dir)

    # Load trees (reconciled)
    T = load_trees(sv.output_dir, prefix=sv.tree_prefix)

    # Test _check_duplication. tiny-phylo uses 'name' not 'recip_paralog'
    p_column = "name"
    expect, obs, df = _check_duplication(sv, T, p_column)
    assert isinstance(expect, int)
    assert isinstance(obs, int)
    assert isinstance(df, pd.DataFrame)

    # Test create_duplications_card
    html = create_duplications_card(sv, T, p_column)
    assert isinstance(html, str)

    # Force a case with excess duplications by creating a mock tree with many 'D' events
    T_copy = T.copy()
    for n in T_copy.traverse():
        if not n.is_leaf:
            n.add_prop("event", "D")
    
    html_warning = create_duplications_card(sv, T_copy, p_column)
    assert "Warning" in html_warning
    assert "We expected" in html_warning