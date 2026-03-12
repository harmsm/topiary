
import pytest
import pandas as pd
import numpy as np
from topiary.ncbi.blast.merge import merge_blast_df, _check_merge, merge_and_annotate

def test__check_merge():
    # 1. merge_list[merge_into] is None
    ml = [None, None, None]
    _check_merge(0, 1, ml)
    assert ml[1] == 1
    assert ml[0] == 1
    
    # 2. merge_list[merge_into] is not None
    ml = [None, 1, None, None]
    _check_merge(2, 1, ml)
    assert ml[2] == 1
    
    # 3. curr_value != merge_into
    ml = [1, 1, 3, 3]
    _check_merge(0, 2, ml)
    assert ml == [3, 3, 3, 3]

def test_merge_blast_df():
    
    # 1. Empty list (hits line 115 after fix)
    result = merge_blast_df([])
    assert result.empty
    
    # 2. Single entry (114)
    df = pd.DataFrame({"accession":["A1"], "query":["Q1"], "subject_start":[1], "subject_end":[10], "e_value":[0], "bits":[100]})
    result = merge_blast_df([df])
    assert len(result) == 1
    
    # 3. Non-overlapping hits (hits line 160)
    data = {"accession":["A1","A1"], "query":["Q1","Q1"], "subject_start":[1, 20], "subject_end":[10, 30], "e_value":[0,0], "bits":[100,100]}
    df_no_overlap = pd.DataFrame(data)
    result = merge_blast_df([df_no_overlap])
    assert len(result) == 2
    
    # 4. Overlapping hits
    data = {
        "accession": ["A1", "A1"],
        "query": ["Q1", "Q2"],
        "subject_start": [1, 5],
        "subject_end": [10, 15],
        "query_start": [1, 1],
        "query_end": [10, 11],
        "bits": [100, 150],
        "e_value": [1e-10, 1e-15],
        "sequence": ["A"*10, "A"*11],
        "hit_def": ["def1", "def1"],
        "hit_id": ["id1", "id1"],
        "title": ["t1", "t1"],
        "length": [100, 100]
    }
    df = pd.DataFrame(data)
    merged = merge_blast_df([df])
    assert len(merged) == 1

def test_merge_and_annotate(mocker):
    
    mocker.patch("topiary.ncbi.blast.merge.get_sequences", return_value=[(None, "M"*100)])
    mocker.patch("topiary.ncbi.blast.merge.parse_ncbi_line", return_value={"iso":1})
    mocker.patch("os.remove")

    # Data with a null title (hits line 241-242)
    data = {
        "accession": ["A1"],
        "query": ["Q1"],
        "subject_start": [1],
        "subject_end": [10],
        "query_start": [1],
        "query_end": [10],
        "bits": [100],
        "e_value": [0],
        "sequence": ["A"*10],
        "hit_def": ["def1"],
        "hit_id": ["id1"],
        "title": [np.nan],
        "length": [100]
    }
    df = pd.DataFrame(data)
    
    out = merge_and_annotate([df], blast_source_list=["local"])
    assert len(out) == 1
