import pytest
from topiary.reports.cards.asr_tree import create_asr_tree_card
from topiary._private.supervisor import Supervisor
from topiary.io.tree import load_trees
import os
import shutil

def test_create_asr_tree_card(tiny_phylo, tmpdir, mocker):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["03_ancestors"]
    sv = Supervisor(calc_dir)

    # Load trees (reconciled)
    T = load_trees(sv.output_dir, prefix=sv.tree_prefix)

    # Output directory for card
    out_dir = str(tmpdir.mkdir("output"))

    # Mock canvas and drawing
    mock_canvas = mocker.Mock()
    mocker.patch("topiary.draw.tree", return_value=mock_canvas)
    
    # Patch canvas_to_html in the asr_tree module
    mocker.patch("topiary.reports.cards.asr_tree.canvas_to_html", return_value="<div>Mock Tree Canvas</div>")

    # Test create_asr_tree_card with ancestor_directory
    html = create_asr_tree_card(sv, out_dir, calc_dir, T)
    assert "Tree used for ASR" in html
    assert "Mock Tree Canvas" in html
    assert os.path.exists(os.path.join(out_dir, "asr-trees.newick"))
    
    # Test with recip_paralog
    sv_recip = Supervisor(calc_dir)
    sv_recip.df["recip_paralog"] = "paralog"
    out_dir_2 = str(tmpdir.mkdir("output2"))
    html = create_asr_tree_card(sv_recip, out_dir_2, None, T)
    assert "Tree used for ASR" in html

    # Test with nickname
    sv_nick = Supervisor(calc_dir)
    # Remove recip_paralog if it was copied (it shouldn't be in sv)
    if "recip_paralog" in sv_nick.df.columns:
        sv_nick.df = sv_nick.df.drop(columns=["recip_paralog"])
    sv_nick.df["nickname"] = "nick"
    out_dir_3 = str(tmpdir.mkdir("output3"))
    html = create_asr_tree_card(sv_nick, out_dir_3, None, T)
    assert "Tree used for ASR" in html

    # Test with default (already covered by first call, but being explicit)
    sv_default = Supervisor(calc_dir)
    if "recip_paralog" in sv_default.df.columns:
        sv_default.df = sv_default.df.drop(columns=["recip_paralog"])
    if "nickname" in sv_default.df.columns:
        sv_default.df = sv_default.df.drop(columns=["nickname"])
    out_dir_4 = str(tmpdir.mkdir("output4"))
    html = create_asr_tree_card(sv_default, out_dir_4, None, T)
    assert "Tree used for ASR" in html