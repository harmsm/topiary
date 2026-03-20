import pytest
from topiary.reports.cards.species_tree import create_species_tree_card
from topiary._private.supervisor import Supervisor
import os

def test_create_species_tree_card(tiny_phylo, tmpdir, mocker):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["03_ancestors"]
    sv = Supervisor(calc_dir)

    # Output directory for card
    out_dir = str(tmpdir.mkdir("output"))

    # Mock canvas and drawing
    mock_canvas = mocker.Mock()
    mocker.patch("topiary.draw.species_tree", return_value=mock_canvas)
    
    # Patch canvas_to_html in the species_tree module
    mocker.patch("topiary.reports.cards.species_tree.canvas_to_html", return_value="<div>Mock Canvas</div>")

    # Test create_species_tree_card
    html = create_species_tree_card(sv, out_dir)
    assert "Species tree for gene/species tree reconciliation" in html
    assert "Mock Canvas" in html
    
    # Verify species-tree.newick was copied
    assert os.path.exists(os.path.join(out_dir, "species-tree.newick"))
    
    import topiary
    topiary.draw.species_tree.assert_called_with(mocker.ANY, output_file=os.path.join(out_dir, "species-tree.pdf"), return_canvas=True)