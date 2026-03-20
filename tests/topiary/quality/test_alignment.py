import pytest
import topiary
from topiary.quality.alignment import score_alignment
import numpy as np
import pandas as pd

def test_score_alignment(test_dataframes, mocker):
    
    # Use a fresh copy for each test case to avoid side effects
    df = test_dataframes["good-df"].copy()
    
    # Test ValueError for invalid align_trim
    with pytest.raises(ValueError):
        score_alignment(df.copy(), align_trim=(1, 0))
    
    with pytest.raises(ValueError):
        score_alignment(df.copy(), align_trim=(0, 1, 2))

    # Test with no alignment column (should trigger topiary.align)
    df_no_align = df.copy()
    if "alignment" in df_no_align.columns:
        df_no_align = df_no_align.drop(columns=["alignment"])
    
    mocker.patch("topiary.align", side_effect=lambda x, **kwargs: x.assign(alignment="AAAAAAA"))
    out_df = score_alignment(df_no_align)
    assert "alignment" in out_df.columns
    
    # Test silent mode
    if "alignment" not in df.columns:
        df["alignment"] = ["AAAAAAA"] * len(df)
        
    out_df = score_alignment(df.copy(), silent=True)
    
    # Test with all gaps column
    df_gaps = df.copy()
    df_gaps["alignment"] = ["---"] * len(df)
    out_df = score_alignment(df_gaps)
    assert "fx_in_sparse" in out_df.columns
    
    # Test with no sparse columns
    out_df = score_alignment(df.copy(), sparse_column_cutoff=1.0)
    assert np.all(out_df.fx_in_sparse.dropna() == 0)
    
    # Test with all sparse columns
    out_df = score_alignment(df.copy(), sparse_column_cutoff=0.0)
    assert "fx_in_sparse" in out_df.columns

    # Test align_trim index forcing (front_index == back_index)
    score_alignment(df.copy(), align_trim=(0.0, 0.000001))
    score_alignment(df.copy(), align_trim=(0.9999999, 1.0))
