
import pytest
import urllib
import http
import io
import pandas as pd
from Bio import Entrez, SeqIO
from topiary.ncbi.entrez.sequences import _get_sequences_thread_function

def test__get_sequences_thread_function_logic(mocker):
    
    import topiary.ncbi.entrez.sequences
    topiary.ncbi.entrez.sequences.topiary.ncbi.NCBI_REQUEST_FREQ = 0
    
    mock_efetch = mocker.patch("Bio.Entrez.efetch")
    mock_handle = mocker.MagicMock()
    # Correct fasta format
    mock_handle.read.return_value = ">ACC1\nMMMMMMMMMM\n"
    mock_efetch.return_value = mock_handle
    
    lock = mocker.MagicMock()
    
    # 1. Success
    out = _get_sequences_thread_function("ACC1", num_tries_allowed=3, lock=lock)
    assert out[0] == ("ACC1", "MMMMMMMMMM")
    
    # 2. HTTPError retry then success (line 55-56, 60-61)
    mock_efetch.side_effect = [urllib.error.HTTPError("url", 404, "Not Found", {}, None), mock_handle]
    out = _get_sequences_thread_function("ACC1", num_tries_allowed=3, lock=lock)
    assert out[0] == ("ACC1", "MMMMMMMMMM")
    
    # 3. Total failure (line 67-68)
    mock_efetch.side_effect = urllib.error.HTTPError("url", 404, "Not Found", {}, None)
    with pytest.raises(ValueError):
        _get_sequences_thread_function("ACC1", num_tries_allowed=1, lock=lock)

    # 4. ID mapping logic (line 86-112)
    # ids = "ACC1,ACC2". returns only ACC1.
    mock_handle.read.return_value = ">ACC1\nMMMMMMMMMM\n"
    mock_efetch.side_effect = [mock_handle]
    out = _get_sequences_thread_function("ACC1,ACC2", num_tries_allowed=1, lock=lock)
    # out should be [("ACC1", "MMMMMMMMMM"), ("ACC2", None)]
    assert out[0] == ("ACC1", "MMMMMMMMMM")
    assert out[1] == ("ACC2", None)
    
    # Complex ID matching (line 90-95)
    # MUST have mismatched counts to hit line 93, 95
    # pdb|1STN|A --> 1STN_A (len 3)
    mock_handle.read.return_value = ">pdb|1STN|A\nMMMMMMMMMM\n"
    mock_efetch.side_effect = [mock_handle]
    out = _get_sequences_thread_function("1STN_A,OTHER", num_tries_allowed=1, lock=lock)
    assert out[0] == ("pdb|1STN|A", "MMMMMMMMMM")
    assert out[1] == ("OTHER", None) # Matches line 95

    # XX|blah.1 --> blah (len 2)
    mock_handle.read.return_value = ">XX|blah.1\nMMMMMMMMMM\n"
    mock_efetch.side_effect = [mock_handle]
    out = _get_sequences_thread_function("blah,OTHER", num_tries_allowed=1, lock=lock)
    assert out[0] == ("XX|blah.1", "MMMMMMMMMM")
    assert out[1] == ("OTHER", None) # Matches line 93
    
    # Trigger ValueError for len > 3 (line 97-98)
    # Need number of ids to be DIFFERENT from number of sequences
    mock_handle.read.return_value = ">A|B|C|D\nMMMMMMMMMM\n"
    mock_efetch.side_effect = [mock_handle]
    with pytest.raises(ValueError):
        # 2 ids, 1 record with len 4 split_list
        _get_sequences_thread_function("X,Y", num_tries_allowed=1, lock=lock)
