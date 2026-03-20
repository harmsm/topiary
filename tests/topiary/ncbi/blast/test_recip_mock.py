
import pytest
import pandas as pd
import numpy as np
import re
from topiary.ncbi.blast.recip import recip_blast, _run_blast

def test_recip_blast_validation():
    
    # 1. Invalid dataframe
    with pytest.raises(ValueError):
        recip_blast("not_a_df", {"P1":["pat1"]})
        
    # 2. No database
    df = pd.DataFrame({
        "accession":["A1"], 
        "sequence":["M"*10], 
        "species":["S1"],
        "name":["N1"],
        "keep":[True]
    })
    with pytest.raises(ValueError):
        recip_blast(df, {"P1":["pat1"]}, local_blast_db=None, ncbi_blast_db=None)

    # 3. TypeError in slicing (line 157-159)
    # Trigger by non-integer start/end
    df_slice = df.copy()
    df_slice["start"] = "not_an_int"
    df_slice["end"] = 10
    with pytest.raises(ValueError) as excinfo:
        recip_blast(df_slice, {"P1":["pat1"]}, local_blast_db="db", use_start_end=True)
    assert "start and end columns must be integers" in str(excinfo.value)

def test__run_blast(mocker):
    # Tests internal _run_blast logic (157-187)
    
    mocker.patch("topiary.ncbi.blast.recip.local_blast", return_value=["local"])
    mocker.patch("topiary.ncbi.blast.recip.ncbi_blast", return_value=["ncbi"])
    
    # 1. Local
    out = _run_blast(["SEQ"], "db", None, None, 10, 0.01, (11,1), 1, False)
    assert out == ["local"]
    
    # 2. NCBI
    out = _run_blast(["SEQ"], None, "nr", None, 10, 0.01, (11,1), 1, False)
    assert out == ["ncbi"]

    # 3. TaxID conflict (line 221-223)
    with pytest.raises(ValueError) as excinfo:
        _run_blast(["SEQ"], None, "nr", 9606, 10, 0.01, (11,1), 1, False, taxid=9606)
    assert "please specify the taxid to use" in str(excinfo.value)

def test_recip_blast_logic(mocker):
    
    # Mock _run_blast
    hit_df = pd.DataFrame({
        "query": ["count0"],
        "accession": ["REC1"],
        "title": ["reciprocal hit to P1"],
        "hit_def": ["reciprocal hit to P1"],
        "bits": [100.0],
        "e_value": [1e-10],
        "subject_start": [1],
        "subject_end": [10]
    })
    mocker.patch("topiary.ncbi.blast.recip._run_blast", return_value=[hit_df])
    
    df = pd.DataFrame({
        "accession": ["A1"],
        "sequence": ["M"*10],
        "species": ["S1"],
        "name": ["A1_name"],
        "keep": [True],
        "uid": ["abcdefghij"]
    })
    
    paralog_patterns = {"P1": ["P1"]}
    
    # 1. Success match
    out = recip_blast(df, paralog_patterns, local_blast_db="db")
    assert out["recip_found_paralog"].iloc[0] == True
    
    # 2. No match found
    paralog_patterns_no_match = {"P2": ["nothing"]}
    out = recip_blast(df, paralog_patterns_no_match, local_blast_db="db")
    assert out["recip_found_paralog"].iloc[0] == False
    
    # 3. No hits from BLAST at all
    mocker.patch("topiary.ncbi.blast.recip._run_blast", return_value=[pd.DataFrame()])
    out = recip_blast(df, paralog_patterns, local_blast_db="db")
    assert out["recip_found_paralog"].iloc[0] == False

    # 4. always_keep renaming (line 486-490)
    # Trigger renaming by having no match (recip_found_paralog is False, recip_paralog is None)
    df_always = df.copy()
    df_always["always_keep"] = True
    df_always["name"] = "test_name"
    # Use a mock that returns nothing to ensure recip_paralog is None
    mocker.patch("topiary.ncbi.blast.recip._run_blast", return_value=[pd.DataFrame()])
    out = recip_blast(df_always, {"P1":["nothing"]}, local_blast_db="db")
    assert out.loc[0, "recip_paralog"] == "test_name"
    assert out.loc[0, "keep"] == True
