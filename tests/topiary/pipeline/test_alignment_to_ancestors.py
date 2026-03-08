import pytest
import topiary
import os
import shutil
import pandas as pd
import numpy as np
from topiary.pipeline.alignment_to_ancestors import _check_restart
from topiary.pipeline.alignment_to_ancestors import alignment_to_ancestors

def test__check_restart(tmpdir, mocker):

    # Case 1: restart = False
    assert _check_restart("some_dir", False) == True

    # Case 2: restart = True, Supervisor fails (file not found)
    mocker.patch("topiary._private.Supervisor", side_effect=FileNotFoundError)
    assert _check_restart("some_dir", True) == True

    # Case 3: restart = True, Supervisor succeeds, calc_status = complete
    mock_sv = mocker.Mock()
    mock_sv.run_parameters = {"calc_status": "complete"}
    mocker.patch("topiary._private.Supervisor", return_value=mock_sv)
    assert _check_restart("some_dir", True) == False

    # Case 4: restart = True, Supervisor succeeds, calc_status != complete
    mock_sv.run_parameters = {"calc_status": "incomplete"}
    mock_rmtree = mocker.patch("topiary.pipeline.alignment_to_ancestors.rmtree")
    mocker.patch("os.path.isdir", return_value=True)
    assert _check_restart("some_dir", True) == True
    mock_rmtree.assert_called_with("some_dir")

def test_alignment_to_ancestors(tmpdir, mocker):
    
    # Mock all the things!
    mock_read_df = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.read_dataframe")
    mock_check_df = mocker.patch("topiary.pipeline.alignment_to_ancestors.check.check_topiary_dataframe")
    mock_check_bool = mocker.patch("topiary.pipeline.alignment_to_ancestors.check.check_bool", side_effect=lambda x, y: x)
    mock_check_int = mocker.patch("topiary.pipeline.alignment_to_ancestors.check.check_int", side_effect=lambda x, y, minimum_allowed=None: x)
    mock_check_float = mocker.patch("topiary.pipeline.alignment_to_ancestors.check.check_float", side_effect=lambda x, y, minimum_allowed=None, maximum_allowed=None: x)
    
    mock_ott_to_mrca = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.opentree.ott_to_mrca", return_value={"is_microbial": False})
    mock_validate_stack = mocker.patch("topiary.pipeline.alignment_to_ancestors.installed.validate_stack")
    mock_check_mpi = mocker.patch("topiary.pipeline.alignment_to_ancestors.check_mpi_configuration")
    mock_df_to_species = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.df_to_species_tree")
    
    mock_find_best_model = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.find_best_model")
    mock_gen_ml_tree = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.generate_ml_tree")
    mock_gen_ancestors = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.generate_ancestors")
    mock_reconcile = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.reconcile")
    mock_gen_bootstraps = mocker.patch("topiary.pipeline.alignment_to_ancestors.topiary.generate_bootstraps")
    mock_pipeline_report = mocker.patch("topiary.pipeline.alignment_to_ancestors.pipeline_report")
    
    mock_check_restart = mocker.patch("topiary.pipeline.alignment_to_ancestors._check_restart", return_value=True)

    # Some setup
    df = pd.DataFrame({"ott": ["ott1", "ott2"], "alignment": ["seq1", "seq2"]})
    mock_check_df.return_value = df
    
    # Create a real tmpdir and change to it
    os.chdir(tmpdir)
    
    # Run simple version (no bootstrap, no reconcile)
    alignment_to_ancestors(df, out_dir="test_out", no_bootstrap=True, force_no_reconcile=True)
    
    assert os.path.isdir("test_out")
    mock_find_best_model.assert_called()
    mock_gen_ml_tree.assert_called()
    mock_gen_ancestors.assert_called()
    assert mock_reconcile.call_count == 0
    assert mock_gen_bootstraps.call_count == 0
    mock_pipeline_report.assert_called()

    # Reset mocks for next run
    mock_find_best_model.reset_mock()
    mock_gen_ml_tree.reset_mock()
    mock_gen_ancestors.reset_mock()
    mock_pipeline_report.reset_mock()

    # Run with reconcile and bootstrap
    alignment_to_ancestors(df, out_dir="test_out_rec", no_bootstrap=False, force_reconcile=True)
    
    assert os.path.isdir("test_out_rec")
    mock_find_best_model.assert_called()
    mock_gen_ml_tree.assert_called()
    assert mock_gen_ancestors.call_count == 2 # Once for ML tree, once for reconciled
    mock_reconcile.assert_called()
    mock_gen_bootstraps.assert_called()
    mock_pipeline_report.assert_called()

    # Test error cases
    from topiary._private.interface import WrappedFunctionException
    
    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, force_reconcile=True, force_no_reconcile=True)
        
    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, horizontal_transfer=True, force_no_reconcile=True)

    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, num_threads=0)

    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, restart=True, overwrite=True)

    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, restart=True, out_dir=None)

    # Test directory exists error
    os.mkdir("exists")
    with pytest.raises(WrappedFunctionException):
        alignment_to_ancestors(df, out_dir="exists", overwrite=False)
