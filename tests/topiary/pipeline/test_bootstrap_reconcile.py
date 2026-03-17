import pytest
import topiary
import os
import glob
from topiary.pipeline.bootstrap_reconcile import bootstrap_reconcile
from topiary._private.interface import WrappedFunctionException

def test_bootstrap_reconcile(tmpdir, mocker):
    
    # Mock all the things!
    mock_check_int = mocker.patch("topiary.pipeline.bootstrap_reconcile.check.check_int", side_effect=lambda x, y, minimum_allowed=None: x)
    mock_check_bool = mocker.patch("topiary.pipeline.bootstrap_reconcile.check.check_bool", side_effect=lambda x, y: x)
    mock_validate_stack = mocker.patch("topiary.pipeline.bootstrap_reconcile.installed.validate_stack")
    mock_check_mpi = mocker.patch("topiary.pipeline.bootstrap_reconcile.check_mpi_configuration")
    
    mock_supervisor_class = mocker.patch("topiary.pipeline.bootstrap_reconcile.Supervisor")
    mock_reconcile = mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary.reconcile")
    mock_reconcile_bootstrap = mocker.patch("topiary.pipeline.bootstrap_reconcile.reconcile_bootstrap")
    mock_pipeline_report = mocker.patch("topiary.pipeline.bootstrap_reconcile.pipeline_report")
    
    # Mock mpi functions
    mock_get_num_slots = mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_num_slots", return_value=56)
    mock_get_hosts = mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_hosts", return_value=["n1"]*56)
    
    # Mock os/glob functions
    mocker.patch("os.path.isdir", return_value=True)
    mocker.patch("glob.glob", side_effect=[
        ["04_bootstraps"], # bootstrap_dirs
        ["rep1.phy", "rep2.phy"] # num_replicates
    ])
    
    # Setup Supervisor mock
    mock_supervisor = mocker.Mock()
    mock_supervisor.status = "complete"
    mock_supervisor.previous_entries = []
    mock_supervisor.run_parameters = {"allow_horizontal_transfer": False}
    mock_supervisor.df = mocker.Mock()
    mock_supervisor.model = "LG+G8"
    mock_supervisor.gene_tree = "((a,b),c);"
    mock_supervisor.species_tree = "((a,b),c);"
    mock_supervisor.reconciled_tree = "((a,b),c);"
    mock_supervisor.seed = 12345
    mock_supervisor_class.return_value = mock_supervisor

    # Define a flexible isdir mock
    isdir_returns = {"prev_dir": True, "05_reconciled-tree-bootstraps": False}
    def side_effect_isdir(path):
        return isdir_returns.get(path, False)
    
    mocker.patch("os.path.isdir", side_effect=side_effect_isdir)
    mocker.patch("os.chdir")

    # Case 1: Normal run (not restart)
    bootstrap_reconcile("prev_dir", num_threads=2)
    
    mock_reconcile.assert_called_with(
        prev_calculation="04_bootstraps",
        calc_dir="05_reconciled-tree-bootstraps",
        bootstrap=True,
        overwrite=False,
        num_threads=2,
        threads_per_rep=8, # Closest factor of 56 to 10
        raxml_binary=mocker.ANY,
        generax_binary=mocker.ANY
    )
    mock_pipeline_report.assert_called()

    # Case 2: Restart
    mock_reconcile.reset_mock()
    isdir_returns["05_reconciled-tree-bootstraps"] = True
    isdir_returns["05_reconciled-tree-bootstraps/working/replicates"] = True
    mocker.patch("glob.glob", side_effect=[
        ["04_bootstraps"], # bootstrap_dirs
        ["rep1.phy", "rep2.phy"] # num_replicates
    ])
    
    bootstrap_reconcile("prev_dir", num_threads=2, restart=True)
    
    mock_reconcile_bootstrap.assert_called()
    assert mock_reconcile.call_count == 0

    # Case 3: Error cases
    
    # previous_run_dir doesn't exist
    isdir_returns["prev_dir"] = False
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile("prev_dir", num_threads=2)

    # No bootstrap dirs
    isdir_returns["prev_dir"] = True
    mocker.patch("glob.glob", return_value=[])
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile("prev_dir", num_threads=2)

    # Supervisor not complete
    mocker.patch("glob.glob", return_value=["04_bootstraps"])
    mock_supervisor.status = "incomplete"
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile("prev_dir", num_threads=2)

    # num_threads > replicates
    mock_supervisor.status = "complete"
    mocker.patch("glob.glob", side_effect=[
        ["04_bootstraps"], # bootstrap_dirs
        ["rep1.phy"] # num_replicates
    ])
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile("prev_dir", num_threads=2)

    # Test threads_per_replicate passed correctly
    isdir_returns["05_reconciled-tree-bootstraps"] = False
    mocker.patch("glob.glob", side_effect=[
        ["04_bootstraps"], # bootstrap_dirs
        ["rep1.phy", "rep2.phy"] # num_replicates
    ])
    mock_reconcile.reset_mock()
    bootstrap_reconcile("prev_dir", num_threads=4, threads_per_replicate=2)
    mock_reconcile.assert_called()
    
    # Test auto-detection
    mock_get_num_slots = mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_num_slots", return_value=56)
    mock_get_hosts = mocker.patch("topiary.pipeline.bootstrap_reconcile.topiary._private.mpi.get_hosts", return_value=["n1"]*56)
    
    mocker.patch("glob.glob", side_effect=[
        ["04_bootstraps"], # bootstrap_dirs
        ["rep"]*100 # many replicates
    ])
    mock_reconcile.reset_mock()
    bootstrap_reconcile("prev_dir") # num_threads=-1, threads_per_replicate=None
    
    # 56 cores -> factor of 56 closest to 10 is 8.
    mock_reconcile.assert_called_with(
        prev_calculation="04_bootstraps",
        calc_dir="05_reconciled-tree-bootstraps",
        bootstrap=True,
        overwrite=False,
        num_threads=56,
        threads_per_rep=8,
        raxml_binary=mocker.ANY,
        generax_binary=mocker.ANY
    )
