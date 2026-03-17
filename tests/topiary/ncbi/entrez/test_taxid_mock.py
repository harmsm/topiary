
import pytest
from topiary.ncbi.entrez.taxid import get_taxid

def test_get_taxid_mock(mocker):
    
    mock_esearch = mocker.patch("Bio.Entrez.esearch")
    mock_read = mocker.patch("Bio.Entrez.read")
    
    # 1. Successful singleton
    mock_read.return_value = {
        "Count": "1",
        "IdList": ["9606"],
        "ErrorList": {}
    }
    result = get_taxid("Homo sapiens")
    assert result == "9606"
    
    # 2. Successful list
    mock_read.return_value = {
        "Count": "2",
        "IdList": ["9606", "10090"],
        "ErrorList": {}
    }
    result = get_taxid(["Homo sapiens", "Mus musculus"])
    assert len(result) == 2
    assert "9606" in result
    assert "10090" in result
    
    # 3. Empty list
    result = get_taxid([])
    assert result == []
    
    # 4. RuntimeError on mismatch
    mock_read.return_value = {
        "Count": "0",
        "IdList": [],
        "ErrorList": {"Field1": ["Error1"]}
    }
    with pytest.raises(RuntimeError) as excinfo:
        get_taxid(["Not a species"])
    assert "Error1" in str(excinfo.value)
    
    # 5. ValueError on bad input types (handled by check.check_iter)
    bad_species = [None, 1.5, {1: "test"}]
    for b in bad_species:
        with pytest.raises(ValueError):
            get_taxid(b)

