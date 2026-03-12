import pytest
from topiary.reports.cards.param import create_param_card
from topiary._private.supervisor import Supervisor
import os

def test_create_param_card(tiny_phylo):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["03_ancestors"]
    sv = Supervisor(calc_dir)

    # anc_dict for testing
    anc_dict = {"anc1": {"bs_support": 100}}

    # Test with ancestor_directory
    html = create_param_card(sv, anc_dict, calc_dir)
    assert "Run parameters" in html
    assert "Evolutionary model" in html
    assert "LG" in html
    assert "Ancestor alt-all cutoff" in html
    assert "0.25" in html

    # Test without ancestor_directory
    html = create_param_card(sv, anc_dict, None)
    assert "Run parameters" in html
    # Check that 0.25 (the cutoff value) is NOT in the HTML when ancestor_directory is None
    assert "0.25" not in html

    # Test with gene tree (not reconciled)
    sv_gene = Supervisor(tiny_phylo["01_gene-tree"])
    html = create_param_card(sv_gene, anc_dict, None)
    assert "Reconciled gene and species tree" in html
    
    # Test with no supports
    anc_dict_no_bs = {"anc1": {"bs_support": None}}
    html = create_param_card(sv, anc_dict_no_bs, None)
    assert "Supports calculated" in html