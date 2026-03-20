import pytest
import topiary
import os
import shutil
import pandas as pd
from topiary.pipeline.seed_to_alignment import _check_restart, _parse_arguments, seed_to_alignment
from topiary._private.interface import WrappedFunctionException

def test__check_restart(tmpdir):
    
    # Case 1: restart = False
    assert _check_restart("some_file", False) == True
    
    # Case 2: restart = True, file exists
    f = tmpdir.join("exists.json")
    f.write("test")
    assert _check_restart(str(f), True) == False
    
    # Case 3: restart = True, file doesn't exist
    assert _check_restart("non_existent.json", True) == True

def test__parse_arguments(mocker):
    
    mock_validate = mocker.patch("topiary.pipeline.seed_to_alignment.installed.validate_stack")
    mock_check_bool = mocker.patch("topiary.pipeline.seed_to_alignment.check.check_bool", side_effect=lambda x, y: x)
    mocker.patch("os.path.exists", return_value=False)
    
    # Normal case
    overwrite, restart, out_dir, species_aware = _parse_arguments(out_dir="test_out")
    assert out_dir == "test_out"
    assert not overwrite
    assert not restart
    
    # Error: overwrite and restart
    with pytest.raises(ValueError):
        _parse_arguments(overwrite=True, restart=True)
        
    # Error: restart without out_dir
    with pytest.raises(ValueError):
        _parse_arguments(restart=True, out_dir=None)

def test_seed_to_alignment(tmpdir, mocker):
    
    # Mock all the things!
    mock_validate = mocker.patch("topiary.pipeline.seed_to_alignment.installed.validate_stack")
    mock_check_bool = mocker.patch("topiary.pipeline.seed_to_alignment.check.check_bool", side_effect=lambda x, y: x)
    mock_read_seed = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.io.read_seed", return_value=(mocker.Mock(), ["spec1"], ["para*"], True))
    mock_df_from_seed = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.df_from_seed", return_value=(pd.DataFrame(), ["spec1"], ["para*"], True))
    mock_write_df = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_dataframe")
    mock_read_df = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.read_dataframe", return_value=pd.DataFrame())
    mock_recip_blast = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.recip_blast", return_value=pd.DataFrame())
    mock_shrink = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.shrink_dataset", return_value=pd.DataFrame())
    mock_align = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.muscle.align", return_value=pd.DataFrame())
    mock_polish = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.polish_alignment", return_value=pd.DataFrame())
    mock_write_fasta = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_fasta")
    
    # Mock NCBI calls
    mock_get_taxid = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_taxid", return_value="12345")
    mock_get_proteome = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_proteome", return_value="proteome.fasta")
    mock_make_blast_db = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.make_blast_db")
    
    # Mock os/shutil
    exists_returns = {"seed.csv": True, "test_out": False}
    def side_effect_exists(path):
        return exists_returns.get(path, False)
    mocker.patch("os.path.exists", side_effect=side_effect_exists)
    
    isdir_returns = {"test_out": False}
    def side_effect_isdir(path):
        return isdir_returns.get(path, False)
    mocker.patch("os.path.isdir", side_effect=side_effect_isdir)
    
    mocker.patch("os.mkdir")
    mocker.patch("os.chdir")
    mocker.patch("shutil.copy")
    
    # Case 1: Normal run from seed_df (str)
    seed_to_alignment("seed.csv", out_dir="test_out")
    
    mock_df_from_seed.assert_called()
    mock_recip_blast.assert_called()
    mock_shrink.assert_called()
    mock_align.assert_called()
    mock_polish.assert_called()
    mock_write_fasta.assert_called()

    # Case 2: Restart (mocking _check_restart return)
    mocker.patch("topiary.pipeline.seed_to_alignment._check_restart", return_value=False)
    seed_to_alignment("seed.csv", out_dir="test_out", restart=True)
    
    assert mock_df_from_seed.call_count == 1 # Still 1 from Case 1
    assert mock_read_df.call_count > 0

    # Error cases (handled by _parse_arguments which we tested above, but some are in seed_to_alignment)
    mocker.patch("topiary.pipeline.seed_to_alignment._check_restart", return_value=True)
    exists_returns["test_out"] = True
    isdir_returns["test_out"] = True
    
    with pytest.raises(WrappedFunctionException):
        seed_to_alignment("seed.csv", out_dir="test_out", overwrite=False, restart=False)
