import pytest
import topiary
import os
import shutil
import json
import pathlib
from pathlib import Path
from topiary.pipeline.bootstrap_reconcile import bootstrap_reconcile
from topiary._private.interface import WrappedFunctionException

def test_restart_logic(tmpdir, mocker):
    """
    Test the restart logic in bootstrap_reconcile when a previous run crashed.
    """
    
    # Mock time-consuming or system-dependent parts
    mocker.patch("topiary.pipeline.bootstrap_reconcile.installed.validate_stack")
    mocker.patch("topiary.pipeline.bootstrap_reconcile.check_mpi_configuration")
    mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_num_slots", return_value=1)
    mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_hosts", return_value=["localhost"])
    mocker.patch("topiary.pipeline.bootstrap_reconcile.pipeline_report")
    
    # Mock the actual core calculation to avoid running generax/raxml
    mock_reconcile_bootstrap = mocker.patch("topiary.pipeline.bootstrap_reconcile.reconcile_bootstrap")
    
    # Create a mock pipeline directory
    pipe_dir = os.path.join(tmpdir, "pipeline")
    os.makedirs(pipe_dir)
    os.chdir(pipe_dir)
    
    # 1. Create a "completed" input bootstrap directory (e.g. from RAxML)
    input_bs_dir = os.path.join(pipe_dir, "04_gene-tree-bootstraps")
    os.makedirs(os.path.join(input_bs_dir, "input"))
    os.makedirs(os.path.join(input_bs_dir, "output", "bootstrap_replicates"))
    
    # Create dummy replicate files so num_replicates > 0
    for i in range(5):
        with open(os.path.join(input_bs_dir, "output", "bootstrap_replicates", f"rep_{i}.phy"), "w") as f:
            f.write("dummy")
            
    # Create run_parameters.json for input_bs_dir with status "complete"
    params = {"calc_status": "complete", "version": topiary.__version__, "seed": 12345}
    with open(os.path.join(input_bs_dir, "run_parameters.json"), "w") as f:
        json.dump(params, f)
        
    # 2. Create a "crashed" reconciliation bootstrap directory
    recon_bs_dir = os.path.join(pipe_dir, "05_reconciled-tree-bootstraps")
    os.makedirs(os.path.join(recon_bs_dir, "working", "replicates", "00001"))
    
    # Create stalled "running" file
    with open(os.path.join(recon_bs_dir, "working", "replicates", "00001", "running"), "w") as f:
        f.write("running")
        
    # Create run_parameters.json for recon_bs_dir with status "running" (crashed)
    params = {
        "calc_status": "running", 
        "version": topiary.__version__, 
        "seed": 12345,
        "calc_type": "reconcile_bootstrap",
        "model": "LG+G8",
        "gene_tree": "input/gene-tree.newick",
        "species_tree": "input/species-tree.newick",
        "reconciled_tree": "input/reconciled-tree.newick",
        "run_parameters": {"allow_horizontal_transfer": True}
    }
    with open(os.path.join(recon_bs_dir, "run_parameters.json"), "w") as f:
        json.dump(params, f)
        
    # Mock Supervisor class to return our dummy attributes
    mock_supervisor_instance = mocker.Mock()
    mock_supervisor_instance.status = "complete" # Default status
    mock_supervisor_instance.previous_entries = []
    mock_supervisor_instance.run_parameters = {"allow_horizontal_transfer": True}
    mock_supervisor_instance.df = mocker.Mock()
    mock_supervisor_instance.model = "LG+G8"
    mock_supervisor_instance.gene_tree = "input/gene-tree.newick"
    mock_supervisor_instance.species_tree = "input/species-tree.newick"
    mock_supervisor_instance.reconciled_tree = "input/reconciled-tree.newick"
    mock_supervisor_instance.seed = 12345
    mock_supervisor_instance.calc_dir = recon_bs_dir # For Case 2
    
    # Define a side effect for Supervisor constructor to return different status
    def supervisor_side_effect(calc_dir=None, seed=None):
        if calc_dir and "04_gene-tree-bootstraps" in calc_dir:
            mock_supervisor_instance.status = "complete"
        elif calc_dir and "05_reconciled-tree-bootstraps" in calc_dir:
            mock_supervisor_instance.status = "running"
        return mock_supervisor_instance

    mocker.patch("topiary.pipeline.bootstrap_reconcile.Supervisor", side_effect=supervisor_side_effect)

    # Now attempt to restart
    # We need to mock os.path.isfile for trees since they don't actually exist
    mocker.patch("os.path.isfile", side_effect=lambda x: True if ".newick" in x or ".phy" in x or "json" in x or "running" in x else os.path.exists(x))

    # Test restart
    bootstrap_reconcile(pipe_dir, restart=True)
    
    # Verify that reconcile_bootstrap was called with restart="replicates"
    called_kwargs = mock_reconcile_bootstrap.call_args[1]
    assert called_kwargs["restart"] == "replicates"
    assert "converge_cutoff" in called_kwargs
    assert called_kwargs["converge_cutoff"] == 0.03
    
    # Verify we are in the pipe_dir or its parent (bootstrap_reconcile changes directory)
    # The function ends by os.chdir('..'), so if it started in pipe_dir, it ends in tmpdir
    assert os.path.abspath(os.getcwd()) == os.path.abspath(tmpdir)

def test_clean_replicate_dir(tmpdir):
    """
    Test that _clean_replicate_dir removes running and skipped files.
    """
    from topiary.generax._reconcile_bootstrap import _clean_replicate_dir
    
    rep_dir = os.path.join(tmpdir, "replicates")
    os.makedirs(os.path.join(rep_dir, "00001"))
    os.makedirs(os.path.join(rep_dir, "00002"))
    
    running1 = os.path.join(rep_dir, "00001", "running")
    skipped2 = os.path.join(rep_dir, "00002", "skipped")
    completed1 = os.path.join(rep_dir, "00001", "completed")
    
    open(running1, 'w').close()
    open(skipped2, 'w').close()
    open(completed1, 'w').close()
    
    _clean_replicate_dir(rep_dir)
    
    assert not os.path.exists(running1)
    assert not os.path.exists(skipped2)
    assert os.path.exists(completed1)
