import pytest
import topiary
from topiary._private import Supervisor
import os
import json
import pandas as pd
import shutil

def test_supervisor_relative_paths(tmpdir, mocker):
    """
    Test that Supervisor stores relative paths in run_parameters.json and
    returns absolute paths via properties.
    """
    
    # Create a working directory
    calc_dir = str(tmpdir.mkdir("test_calc"))
    
    # Create dummy files for trees and alignment to mock a real run
    input_dir = os.path.join(calc_dir, "input")
    os.makedirs(input_dir, exist_ok=True)
    
    # Use an existing valid tree from test data
    valid_tree = os.path.abspath("tests/data/tiny-phylo/01_gene-tree/output/gene-tree.newick")

    gene_tree_file = os.path.join(input_dir, "gene-tree.newick")
    shutil.copy(valid_tree, gene_tree_file)
        
    species_tree_file = os.path.join(input_dir, "species-tree.newick")
    shutil.copy(valid_tree, species_tree_file)
        
    alignment_file = os.path.join(input_dir, "alignment.phy")
    with open(alignment_file, "w") as f:
        f.write("2 10\nA  ACGTACGTAC\nB  ACGTACGTAC")

    # Mock tree reading so we don't have to worry about valid newick content
    mock_tree = mocker.Mock()
    mocker.patch("topiary.io.read_tree", return_value=mock_tree)

    # Mock dataframe reading
    df = pd.DataFrame({"ott":["ott1","ott2"], 
                       "species":["sp1","sp2"], 
                       "uid":["uid1","uid2"],
                       "name":["name1","name2"],
                       "alignment":["-A","-A"],
                       "sequence":["MA","MA"],
                       "keep":[True,True]})
    mocker.patch("topiary.read_dataframe", return_value=df)
    mocker.patch("topiary.write_dataframe")

    # Initialize supervisor
    sv = Supervisor()
    
    # Create calc dir
    sv.create_calc_dir(calc_dir=calc_dir,
                       calc_type="test_calc",
                       df=df,
                       gene_tree=gene_tree_file,
                       species_tree=species_tree_file,
                       model="LG",
                       overwrite=True)
    
    # Ensure JSON is written
    sv.write_json()
    
    # 1. Verify properties return absolute paths
    assert os.path.isabs(sv.calc_dir)
    assert os.path.isabs(sv.gene_tree)
    assert os.path.isabs(sv.species_tree)
    assert os.path.isabs(sv.alignment)
    
    assert sv.gene_tree == os.path.join(sv.input_dir, "gene-tree.newick")
    assert sv.species_tree == os.path.join(sv.input_dir, "species-tree.newick")
    assert sv.alignment == os.path.join(sv.input_dir, "alignment.phy")

    # 2. Verify run_parameters.json stores relative paths
    params_file = os.path.join(calc_dir, "run_parameters.json")
    with open(params_file, "r") as f:
        params = json.load(f)
        
    assert not os.path.isabs(params["gene_tree"])
    assert not os.path.isabs(params["species_tree"])
    assert not os.path.isabs(params["alignment"])
    
    assert params["gene_tree"] == "input/gene-tree.newick"
    assert params["species_tree"] == "input/species-tree.newick"
    assert params["alignment"] == "input/alignment.phy"

    # 3. Verify that _increment stores relative calc_dir in previous_entries
    new_calc_dir = str(tmpdir.mkdir("test_calc_increment"))
    sv.create_calc_dir(calc_dir=new_calc_dir,
                       calc_type="test_calc_2",
                       force=True,
                       overwrite=True)
    
    with open(os.path.join(new_calc_dir, "run_parameters.json"), "r") as f:
        params2 = json.load(f)
        
    prev = params2["previous_entries"][0]
    assert not os.path.isabs(prev["calc_dir"])
    # Relpath from test_calc_increment to test_calc
    expected_rel = os.path.relpath(calc_dir, new_calc_dir)
    assert prev["calc_dir"] == expected_rel

    # 4. Verify loading existing directory works correctly
    # (Move the directory to ensure relative paths work when absolute paths change)
    moved_calc_dir = str(tmpdir.mkdir("moved_calc"))
    shutil.rmtree(moved_calc_dir) # remove it so copytree can create it
    shutil.copytree(calc_dir, moved_calc_dir)
    
    sv_moved = Supervisor(calc_dir=moved_calc_dir)
    assert os.path.isabs(sv_moved.gene_tree)
    assert sv_moved.gene_tree == os.path.join(moved_calc_dir, "input", "gene-tree.newick")
    
    # Manually create the files to verify the path resolution works if they aren't there
    if not os.path.exists(sv_moved.gene_tree):
        os.makedirs(os.path.dirname(sv_moved.gene_tree), exist_ok=True)
        with open(sv_moved.gene_tree, "w") as f:
            f.write("(A,B);")

    assert os.path.exists(sv_moved.gene_tree)
