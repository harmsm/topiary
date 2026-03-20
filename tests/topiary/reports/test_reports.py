import pytest
import os
import shutil
import pandas as pd
import numpy as np
import json

import topiary
from topiary.reports.reports import _find_directories
from topiary.reports.reports import tree_report
from topiary.reports.reports import pipeline_report

def test__find_directories(tmpdir, small_phylo):
    cwd = os.getcwd()
    
    # 3. single calculation directory (DO THIS BEFORE CHDIR TO TMPDIR)
    # This ensures absolute paths work correctly.
    out = _find_directories(small_phylo["00_find-model"])
    assert out['model'] is not None

    os.chdir(tmpdir)

    # Complete directory
    os.mkdir("test1")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors",
                    "03_reconciled-tree",
                    "04_reconciled-tree-ancestors",
                    "05_gene-tree-bootstraps",
                    "06_reconciled-tree-bootstraps"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d], os.path.join("test1", d))
    
    out = _find_directories("test1")

    assert out['model'] == os.path.join('test1', '00_find-model')
    assert out['gene']['anc'] == os.path.join('test1', '02_gene-tree-ancestors')
    assert out['gene']['tree'] == os.path.join('test1', '05_gene-tree-bootstraps')
    assert out['reconciled']['anc'] == os.path.join('test1', '04_reconciled-tree-ancestors')
    assert out['reconciled']['tree'] == os.path.join('test1', '06_reconciled-tree-bootstraps')

    # 1. run_parameters.json missing creation_time
    os.mkdir("test_missing_time")
    shutil.copytree(small_phylo["00_find-model"], os.path.join("test_missing_time", "00_find-model"))
    param_path = os.path.join("test_missing_time", "00_find-model", "run_parameters.json")
    with open(param_path, 'r') as f:
        params = json.load(f)
    if "creation_time" in params:
        del params["creation_time"]
    params["calc_status"] = "complete"
    with open(param_path, 'w') as f:
        json.dump(params, f)
    
    out = _find_directories("test_missing_time")
    assert out['model'] == os.path.join('test_missing_time', '00_find-model')

    # 2. calc_status != "complete"
    os.mkdir("test_incomplete")
    shutil.copytree(small_phylo["00_find-model"], os.path.join("test_incomplete", "calc"))
    param_path = os.path.join("test_incomplete", "calc", "run_parameters.json")
    with open(param_path, 'r') as f:
        params = json.load(f)
    params["calc_status"] = "running"
    with open(param_path, 'w') as f:
        json.dump(params, f)
    
    out = _find_directories("test_incomplete")
    assert out['model'] is None

    # 4. tree_prefix is None
    os.mkdir("test_no_prefix")
    shutil.copytree(small_phylo["00_find-model"], os.path.join("test_no_prefix", "calc"))
    param_path = os.path.join("test_no_prefix", "calc", "run_parameters.json")
    with open(param_path, 'r') as f:
        params = json.load(f)
    params["calc_type"] = "not_a_tree_calc"
    params["calc_status"] = "complete"
    with open(param_path, 'w') as f:
        json.dump(params, f)
    out = _find_directories("test_no_prefix")
    assert out['gene']['tree'] is None

    os.chdir(cwd)

def test_tree_report(tmpdir, small_phylo, mocker):
    
    output_dir = os.path.join(tmpdir, "test_report")
    tree_dir = small_phylo["05_gene-tree-bootstraps"]
    anc_dir = small_phylo["04_reconciled-tree-ancestors"]
    
    import ete4 as ete
    T = ete.Tree("((A:0.1,B:0.1)90:0.1,C:0.2);")
    for i, leaf in enumerate(T.leaves()):
        leaf.add_prop("uid", f"uid{i}")
        leaf.name = f"uid{i}"
        leaf.add_prop("ott", "ott1")
    for node in T.traverse():
        if not node.is_leaf:
            node.add_prop("anc_label", "a0")
            node.add_prop("bs_support", 100)
    
    mocker.patch("topiary.reports.reports.load_trees", return_value=T)
    mocker.patch("topiary.reports.reports.create_main_html", return_value=("top", "bottom"))
    mocker.patch("topiary.reports.reports.create_card", return_value="card")
    mocker.patch("topiary.reports.reports.create_row", return_value="row")
    mocker.patch("topiary.reports.reports.create_icon_row", return_value="icons")
    mocker.patch("topiary.reports.reports.create_info_modal", return_value="modal")
    mocker.patch("topiary.reports.reports.create_ancestor_card", return_value="anc_card")
    mocker.patch("topiary.reports.reports.create_input_card", return_value="input_card")
    mocker.patch("topiary.reports.reports.create_model_card", return_value="model_card")
    mocker.patch("topiary.reports.reports.create_duplications_card", return_value="dup_card")
    mocker.patch("topiary.reports.reports.create_species_tree_card", return_value="species_card")
    mocker.patch("topiary.reports.reports.create_asr_tree_card", return_value="asr_card")
    mocker.patch("topiary.reports.reports.create_param_card", return_value="param_card")
    mocker.patch("topiary.opentree.ott_to_mrca", return_value={"ott_name": "Test Group"})

    from topiary._private import Supervisor
    mock_sv = mocker.MagicMock(spec=Supervisor)
    mock_sv.tree_prefix = "reconciled"
    mock_sv.output_dir = tree_dir
    mock_sv.df = pd.DataFrame({"uid":["uid0","uid1","uid2"],"ott":["ott1","ott1","ott1"],"recip_paralog":["P1","P1","P1"]})
    mock_sv.model = "JTT"
    mocker.patch("topiary.reports.reports.Supervisor", return_value=mock_sv)

    # Run for reconciled tree with zip AND ancestor_directory
    tree_report(tree_directory=tree_dir,
                output_directory=output_dir,
                ancestor_directory=anc_dir,
                overwrite=True,
                create_zip_file=True)
    
    assert os.path.exists(os.path.join(output_dir, "index.html"))
    assert os.path.exists(f"{output_dir}.zip")

def test_pipeline_report(tmpdir, small_phylo, mocker):
    
    pipeline_dir = os.path.join(tmpdir, "pipeline")
    os.mkdir(pipeline_dir)
    dirs_to_copy = ["00_find-model", "05_gene-tree-bootstraps", "06_reconciled-tree-bootstraps"]
    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d], os.path.join(pipeline_dir, d))
    
    output_dir = os.path.join(tmpdir, "pipeline_report")
    
    mocker.patch("topiary.reports.reports.create_main_html", return_value=("top", "bottom"))
    mocker.patch("topiary.reports.reports.create_card", return_value="card")
    mocker.patch("topiary.reports.reports.create_row", return_value="row")
    mocker.patch("topiary.reports.reports.create_icon_row", return_value="icons")
    mocker.patch("topiary.reports.reports.create_info_modal", return_value="modal")
    mocker.patch("topiary.reports.reports.create_ancestor_card", return_value="anc_card")
    mocker.patch("topiary.reports.reports.create_input_card", return_value="input_card")
    mocker.patch("topiary.reports.reports.create_model_card", return_value="model_card")
    mocker.patch("topiary.reports.reports.create_duplications_card", return_value="dup_card")
    mocker.patch("topiary.reports.reports.create_species_tree_card", return_value="species_card")
    mocker.patch("topiary.reports.reports.create_asr_tree_card", return_value="asr_card")
    mocker.patch("topiary.reports.reports.create_param_card", return_value="param_card")
    mocker.patch("topiary.opentree.ott_to_mrca", return_value={"ott_name": "Test Group"})

    import ete4 as ete
    T = ete.Tree("((A:0.1,B:0.1)90:0.1,C:0.2);")
    for i, leaf in enumerate(T.leaves()):
        leaf.add_prop("uid", f"uid{i}")
        leaf.name = f"uid{i}"
        leaf.add_prop("ott", "ott1")
    for node in T.traverse():
        if not node.is_leaf:
            node.add_prop("anc_label", "a0")
            node.add_prop("bs_support", 100)
    
    mocker.patch("topiary.reports.reports.load_trees", return_value=T)
    
    from topiary._private import Supervisor
    
    def mock_sv_init(d):
        sv = mocker.MagicMock(spec=Supervisor)
        sv.status = "complete"
        sv.output_dir = d
        if "reconciled" in d:
            sv.tree_prefix = "reconciled"
        elif "gene" in d or "ml" in d:
            sv.tree_prefix = "gene"
        else:
            sv.tree_prefix = None
        sv.df = pd.DataFrame({"uid":["uid0","uid1","uid2"],"ott":["ott1","ott1","ott1"],"name":["N1","N2","N3"]})
        sv.model = "JTT"
        sv.calc_type = "find_best_model" if "00_find-model" in d else "ml_tree"
        return sv

    mocker.patch("topiary.reports.reports.Supervisor", side_effect=mock_sv_init)

    pipeline_report(pipeline_directory=pipeline_dir, output_directory=output_dir, create_zip_file=True)
    
    assert os.path.exists(os.path.join(output_dir, "index.html"))
    assert os.path.exists(f"{output_dir}.zip")

    empty_dir = os.path.join(tmpdir, "empty")
    os.mkdir(empty_dir)
    with pytest.raises(ValueError):
        pipeline_report(pipeline_directory=empty_dir, output_directory=output_dir)