
import pytest
import os
from topiary.ncbi.entrez.download import ncbi_ftp_download, _read_md5_file

def test__read_md5_file(tmpdir):
    # Hits line 30 (empty line)
    md5_path = tmpdir.join("md5.txt")
    md5_path.write("hash1  ./file1\n\nhash2  ./file2\n")
    d = _read_md5_file(str(md5_path))
    assert d["file1"] == "hash1"
    assert d["file2"] == "hash2"

def test_ncbi_ftp_download_logic(mocker):
    
    # Mock ftp_download to avoid actual FTP
    mock_ftp = mocker.patch("topiary.ncbi.entrez.download.ftp_download")
    mocker.patch("os.path.isfile", return_value=True)
    mocker.patch("os.remove")
    
    # Correct filename based on URL:
    url = "ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"
    fname = "GCF_000001405.40_GRCh38.p14_protein.faa.gz"
    
    # Mock _read_md5_file
    mocker.patch("topiary.ncbi.entrez.download._read_md5_file", return_value={fname:"hash1"})
    
    # Mock calc_md5
    mocker.patch("topiary.ncbi.entrez.download.calc_md5", return_value="hash1")
    
    # 1. Success
    ncbi_ftp_download(url)
    
    # 2. md5 download failure (line 85-89)
    mock_ftp.side_effect = Exception("FTP Fail")
    with pytest.raises(Exception):
        ncbi_ftp_download(url)
        
    # 3. md5 mismatch or file not in md5 (line 93-96)
    mock_ftp.side_effect = None
    mocker.patch("topiary.ncbi.entrez.download._read_md5_file", return_value={})
    with pytest.raises(FileNotFoundError):
        ncbi_ftp_download(url)

def test_ncbi_ftp_download_retry(mocker):
    
    # Mock for main loop retry (line 107-108 etc)
    mock_ftp = mocker.patch("topiary.ncbi.entrez.download.ftp_download")
    mocker.patch("os.path.isfile", return_value=True)
    mocker.patch("os.remove")
    
    url = "ftp.ncbi.nlm.nih.gov/all/GCF_001"
    fname = "GCF_001_protein.faa.gz"
    
    mocker.patch("topiary.ncbi.entrez.download._read_md5_file", return_value={fname:"hash1"})
    mocker.patch("topiary.ncbi.entrez.download.calc_md5", side_effect=["bad_hash", "hash1"])
    
    # This should retry because of bad_hash once, then succeed
    ncbi_ftp_download(url, num_attempts=2)
    assert mock_ftp.call_count >= 2

    # Test file not found after download (line 106-108)
    mocker.patch("os.path.isfile", return_value=False)
    with pytest.raises(RuntimeError):
        ncbi_ftp_download(url, num_attempts=1)
