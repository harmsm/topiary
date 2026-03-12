
import pytest
import pandas as pd
import os
from topiary.ncbi.entrez.proteome import get_proteome

def test_get_proteome_logic(mocker):
    
    # Mock get_taxid to avoid taxonomy search.
    mocker.patch("topiary.ncbi.get_taxid", return_value="9606")
    mock_download = mocker.patch("topiary.ncbi.entrez.proteome.ncbi_ftp_download")
    mocker.patch("os.path.isfile", return_value=True)
    mock_remove = mocker.patch("os.remove")
    
    # Bio.Entrez mocks
    mock_esearch = mocker.patch("Bio.Entrez.esearch")
    mock_esummary = mocker.patch("Bio.Entrez.esummary")
    mock_read = mocker.patch("Bio.Entrez.read")
    
    # Mock records for esummary
    mock_record = {
        "LastUpdateDate": "2023-01-01",
        "RefSeq_category": "reference genome",
        "FtpPath_RefSeq": "ftp://path/GCF_001",
        "FtpPath_GenBank": "ftp://path/GCA_001"
    }
    
    esummary_resp = {
        "DocumentSummarySet": {
            "DocumentSummary": [mock_record]
        }
    }
    
    # 1. Success path
    mock_read.side_effect = [
        {"IdList":["123"]}, # get_proteome_ids Loop 1
        {"IdList":[]},      # get_proteome_ids Loop 2
        esummary_resp       # _get_records call
    ]
    val = get_proteome(species="Homo sapiens")
    assert val is not None

    # 2. FileNotFoundError (line 224-245)
    # This should hit the extensive error message and raise RuntimeError
    mock_download.side_effect = FileNotFoundError("Missing")
    mock_read.side_effect = [
        {"IdList":["123"]},
        {"IdList":[]},
        esummary_resp
    ]
    with pytest.raises(RuntimeError) as excinfo:
        # Pass species ONLY to hit line 227
        get_proteome(species="Homo sapiens")
    assert "Could not download proteome" in str(excinfo.value)
    assert "Homo sapiens" in str(excinfo.value)

    # 3. FileNotFoundError with taxid (line 229)
    # Reset mocks
    mock_download.side_effect = FileNotFoundError("Missing")
    mock_read.side_effect = [
        {"IdList":["123"]},
        {"IdList":[]},
        esummary_resp
    ]
    with pytest.raises(RuntimeError) as excinfo:
        # Pass taxid ONLY to hit line 229
        get_proteome(taxid="9606")
    assert "taxid = 9606" in str(excinfo.value)

    # 4. KeyError in get_proteome_ids (line 133-134)
    # Trigger by missing IdList
    mock_download.side_effect = None
    mock_read.side_effect = [
        {}, # No IdList in first loop
        {"IdList":[]}, # Empty in second
        esummary_resp
    ]
    with pytest.raises(RuntimeError):
         get_proteome(species="Homo sapiens")

    # 5. KeyError in _get_records (line 165-167)
    mock_read.side_effect = [
        {"IdList":["123"]},
        {"IdList":[]},
        {} # Missing DocumentSummarySet
    ]
    with pytest.raises(RuntimeError):
        get_proteome(species="Homo sapiens")

    # 6. Total failure (return None at 259)
    # Trigger by ncbi_ftp_download raising RuntimeError (caught/continued at 222)
    mock_download.side_effect = RuntimeError("Other error")
    mock_read.side_effect = [
        {"IdList":["123"]},
        {"IdList":[]},
        esummary_resp
    ]
    val = get_proteome(species="Homo sapiens")
    assert val is None
    
    # 7. os.remove FileNotFoundError (line 253-254)
    mock_download.side_effect = None
    mock_remove.side_effect = FileNotFoundError 
    mock_read.side_effect = [
        {"IdList":["123"]},
        {"IdList":[]},
        esummary_resp
    ]
    val = get_proteome(species="Homo sapiens")
    assert val is not None
