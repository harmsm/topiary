import pytest
import topiary
from topiary.quality.polish import _get_cutoff, polish_alignment
import numpy as np
import pandas as pd

def test__get_cutoff():
    
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    
    # Test default percentile
    assert _get_cutoff(x, pct=0.975) == 10
    
    # Test with avg_bin_contents
    val = _get_cutoff(x, avg_bin_contents=5, pct=0.9)
    assert np.isclose(val, 7.75)
    
    with pytest.raises(ValueError):
        _get_cutoff(x, avg_bin_contents=0)

def test_polish_alignment(test_dataframes, mocker):
    
    df = test_dataframes["good-df"].copy()
    df["always_keep"] = False
    df.loc[0, "always_keep"] = True
    
    scored_df = df.copy()
    scored_df["fx_in_sparse"] = np.linspace(0, 1, len(df))
    scored_df["sparse_run_length"] = np.linspace(0, 10, len(df))
    scored_df["fx_missing_dense"] = np.linspace(0, 1, len(df))
    scored_df["partial"] = [False] * len(df)
    scored_df.loc[len(df)-1, "partial"] = True
    
    mocker.patch("topiary.quality.polish.score_alignment", return_value=scored_df)
    mock_muscle = mocker.patch("topiary.muscle.align", side_effect=lambda x: x)
    
    # Test polish_alignment
    out_df = polish_alignment(df.copy(), realign=True)
    assert len(out_df) == len(df)
    assert mock_muscle.called
    
    # Test with top_sparse_run == 0 edge case
    scored_df_zero = scored_df.copy()
    scored_df_zero["sparse_run_length"] = 0
    scored_df_zero["keep"] = True
    mocker.patch("topiary.quality.polish.score_alignment", return_value=scored_df_zero)
    
    out_df_zero = polish_alignment(df.copy(), realign=False)
    assert out_df_zero.loc[len(df)-1, "keep"] == False
    
    # Test realign=False
    mock_muscle.reset_mock()
    mocker.patch("topiary.quality.polish.score_alignment", return_value=scored_df)
    out_df_no_realign = polish_alignment(df.copy(), realign=False)
    assert not mock_muscle.called
