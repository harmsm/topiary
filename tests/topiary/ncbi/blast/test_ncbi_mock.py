
import pytest
import os
import pandas as pd
import numpy as np
import urllib
import http
import topiary
from topiary.ncbi.blast.ncbi import ncbi_blast, _prepare_for_blast, _construct_args, _ncbi_blast_thread_function

def test__prepare_for_blast():
    
    # db validation
    with pytest.raises(ValueError):
        _prepare_for_blast("SEQ", db=None, taxid=9606, blast_program="blastp", hitlist_size=10, e_value_cutoff=0.01, gapcosts=(11,1), url_base="https://base.com", kwargs={})

    # taxid validation
    _prepare_for_blast("SEQ", db="nr", taxid="9606", blast_program="blastp", hitlist_size=10, e_value_cutoff=0.01, gapcosts=(11,1), url_base="https://base.com", kwargs={})
    _prepare_for_blast("SEQ", db="nr", taxid=[9606], blast_program="blastp", hitlist_size=10, e_value_cutoff=0.01, gapcosts=(11,1), url_base="https://base.com", kwargs={})

def test__construct_args():
    
    # Sequence too long error
    seq_list_long = ["A" * 200]
    with pytest.raises(ValueError):
        _construct_args(seq_list_long, {"db":"nr"}, max_query_length=100, num_tries_allowed=3, keep_blast_xml=False, num_threads=4)

    # Splitting logic (line 288)
    seq_list_mix = ["A"*30, "A"*30, "A"*30]
    # new_sequence (~38). max_query_length=50.
    # Block 1 starts empty. Seq0 (38 < 50) added. current=38.
    # Seq1 (38+38=76 > 50) -> hits 288.
    _construct_args(seq_list_mix, {"db":"nr"}, max_query_length=50, num_tries_allowed=3, keep_blast_xml=False, num_threads=2)

    # Windowing logic coverage
    _construct_args(["S"]*11, {"db":"nr"}, max_query_length=1000, num_tries_allowed=3, keep_blast_xml=False, num_threads=3)
    _construct_args(["S"]*5, {"db":"nr"}, max_query_length=1000, num_tries_allowed=3, keep_blast_xml=False, num_threads=2)

def test__ncbi_blast_thread_function(mocker):
    
    import topiary.ncbi
    topiary.ncbi.NCBI_REQUEST_FREQ = 0
    
    mock_qblast = mocker.patch("Bio.Blast.NCBIWWW.qblast")
    mock_handle = mocker.MagicMock()
    mock_handle.read.return_value = "<xml>fake</xml>"
    mock_qblast.return_value = mock_handle
    
    df = pd.DataFrame({"query":["count0"], "accession":["ACC1"]})
    mocker.patch("topiary.ncbi.blast.ncbi.read_blast_xml", return_value=([df], ["fake.xml"]))
    mocker.patch("os.remove")
    
    lock = mocker.MagicMock()
    this_query = {"database":"nr", "sequence":">count0\nACGT\n"}

    # Success
    out = _ncbi_blast_thread_function(this_query, num_tries_allowed=3, keep_blast_xml=False, lock=lock)
    assert out.equals(df)

    # CPU limit handling
    mocker.patch("topiary.ncbi.blast.ncbi.read_blast_xml", return_value=(None, ["fake.xml"]))
    mock_qblast.side_effect = [mock_handle]
    assert _ncbi_blast_thread_function(this_query, num_tries_allowed=1, keep_blast_xml=False, lock=lock) is None

def test_ncbi_blast(mocker):
    
    mocker.patch("topiary.ncbi.blast.ncbi.threads.get_num_threads", return_value=2)
    df = pd.DataFrame({"query":["count0"], "accession":["ACC1"]})
    
    def side_effect(kwargs_list, func, num_threads, **kwargs):
        if num_threads > 1:
            return [None] * len(kwargs_list)
        else:
            return [df]
            
    mocker.patch("topiary._private.threads.thread_manager", side_effect=side_effect)
    
    # Run and trigger CPU limit fallback
    result = ncbi_blast(["ACGT", "ACGT"], db="nr", taxid=9606, num_threads=2)
    assert isinstance(result, list)

def test_combine_hits_edge_cases():
    from topiary.ncbi.blast.ncbi import _combine_hits
    
    df_no_hits = pd.DataFrame({"query":["count0"], "accession":[np.nan]})
    result = _combine_hits([df_no_hits], return_singleton=True)
    assert result.empty
