import pytest
from topiary.reports.cards.ancestor import create_ancestor_card
from topiary._private.supervisor import Supervisor
import os
import shutil
import pandas as pd
import numpy as np

def test_create_ancestor_card(tiny_phylo, tmpdir, mocker):
    
    # Use actual tiny_phylo dir as ancestor_directory
    ancestor_directory = tiny_phylo["03_ancestors/"]

    # Output directory for card results
    out_dir = str(tmpdir.mkdir("output"))

    # anc_dict for testing based on tiny_phylo data
    anc_dict = {"anc1": {
        "bs_support": 100,
        "taxonomic_dist": "TestTaxa",
        "paralog_call": "P1",
        "event": "S",
        "paralogs": "P1: 100%",
        "anc_pp": 0.95,
        "num_descendants": 10
    }}

    # event_color and event_name
    event_color = {"S":"#023E55","D":"#64007F","L":"#BAD316","T":"#407E98",None:"#000000"}
    event_name = {"S":"speciation","D":"duplication","L":"loss","T":"transfer",None:'N/A'}

    # Mock plot_ancestor_data and plt.close to avoid heavy drawing
    mock_fig = mocker.Mock()
    mock_ax = mocker.Mock()
    mocker.patch("topiary.reports.cards.ancestor.plot_ancestor_data", return_value=(mock_fig, mock_ax))
    mocker.patch("topiary.reports.cards.ancestor.plt.close")

    # Mock shutil.copy to avoid polluting output directory with large files
    mocker.patch("shutil.copy")

    # Test create_ancestor_card
    html = create_ancestor_card(anc_dict, out_dir, ancestor_directory, event_color, event_name)
    assert "Ancestors" in html
    assert "anc1" in html
    assert "TestTaxa" in html

    # Verify anc1.fasta and anc1.csv were created in out_dir
    assert os.path.exists(os.path.join(out_dir, "anc1.fasta"))
    assert os.path.exists(os.path.join(out_dir, "anc1.csv"))

    # Verify savefig was called for svg and pdf
    # anc1.svg and anc1_pp.pdf
    mock_fig.savefig.assert_any_call(os.path.join(out_dir, "anc1.svg"), bbox_inches="tight")
    mock_fig.savefig.assert_any_call(os.path.join(out_dir, "anc1_pp.pdf"), bbox_inches="tight")