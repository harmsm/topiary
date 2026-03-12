import pytest
from topiary.reports.cards.model import create_model_card
from topiary._private.supervisor import Supervisor
import os
import shutil

def test_create_model_card(tiny_phylo, tmpdir):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["00_find-best-model"]
    sv = Supervisor(calc_dir)

    # Output directory for card
    out_dir = str(tmpdir.mkdir("output"))

    # Test create_model_card
    html = create_model_card(sv, out_dir)
    assert "Evolutionary model selection" in html
    assert "Best model: LG" in html
    
    # Verify model-comparison.csv was copied
    assert os.path.exists(os.path.join(out_dir, "model-comparison.csv"))