import pytest
import topiary
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt

from topiary.draw.ancestor_data import _draw_histogram
from topiary.draw.ancestor_data import plot_ancestor_data

def test__draw_histogram():
    fig, ax = plt.subplots()
    values = np.random.random(100)
    max_counts = _draw_histogram(values, ax)
    
    assert max_counts > 0
    assert len(ax.patches) > 0
    plt.close(fig)

def test_plot_ancestor_data():
    # Create dummy ancestral dataframe
    sites = np.arange(10)
    ml_pp = np.random.random(10)
    alt_pp = np.random.random(10) * 0.5
    site_type = ["good"] * 8 + ["gap"] * 1 + ["possible gap"] * 1
    
    df_anc = pd.DataFrame({
        "site": sites,
        "ml_pp": ml_pp,
        "alt_pp": alt_pp,
        "site_type": site_type
    })
    
    # Test plotting with various options
    fig, ax = plot_ancestor_data(df_anc, anc_name="test_anc", anc_data_string="test data")
    
    assert fig is not None
    assert len(ax) == 2
    assert os.path.exists("test_anc.pdf")
    
    # Clean up
    os.remove("test_anc.pdf")
    plt.close(fig)

    # Test with close_plot=True
    fig, ax = plot_ancestor_data(df_anc, close_plot=True)
    # Fig should be closed, but we can still check properties if not fully purged
    assert fig is not None 

    # Test with contiguous gaps
    df_anc.loc[5:7, "site_type"] = "gap"
    fig, ax = plot_ancestor_data(df_anc, close_plot=True)

def test_plot_ancestor_data_no_gaps():
    # sites with no gaps (sites = 0,1,2, nogap = 0,1,2, gap = empty)
    sites = np.arange(3)
    ml_pp = np.array([0.9, 0.9, 0.9])
    alt_pp = np.array([0.1, 0.1, 0.1])
    site_type = ["good", "good", "good"]
    
    df_anc = pd.DataFrame({
        "site": sites,
        "ml_pp": ml_pp,
        "alt_pp": alt_pp,
        "site_type": site_type
    })
    
    fig, ax = plot_ancestor_data(df_anc, close_plot=True)
    assert fig is not None
