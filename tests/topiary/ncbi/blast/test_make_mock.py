
import pytest
import os
import subprocess
import gzip
import io
from topiary.ncbi.blast.make import make_blast_db

def test_make_blast_db(tmpdir, mocker):

    # Preparation
    os.chdir(tmpdir)
    faa1 = tmpdir.join("test1.faa")
    faa1.write(">seq1\nACGT\n")
    
    faa2 = tmpdir.join("test2.faa.gz")
    with gzip.open(str(faa2), "wb") as f:
        f.write(b">seq2\nACGT\n")
    
    bad_ext = tmpdir.join("test3.txt")
    bad_ext.write("not a faa")

    # Mocks
    mock_run = mocker.patch("subprocess.run")
    mock_popen = mocker.patch("subprocess.Popen")
    mock_popen_instance = mock_popen.return_value
    mock_popen_instance.stdout = io.StringIO("building...")
    mock_popen_instance.wait.return_value = 0

    # 1. Basic success with multiple files
    make_blast_db([str(faa1), str(faa2)], "db_name")
    assert mock_run.called
    assert mock_popen.called
    assert os.path.exists("db_name") == False # It just creates the cmd, mock doesn't create files

    # Reset mock stdout for next call
    mock_popen_instance.stdout = io.StringIO("building...")
    # 2. String input instead of list
    make_blast_db(str(faa1), "db_name_str")
    
    # 3. Invalid input_files type
    with pytest.raises(ValueError):
        make_blast_db(123, "db")

    # 4. Overwrite behavior
    # Manually create ONLY .psq to trigger FileNotFoundError on others
    with open("db_exists.psq", "w") as f:
        f.write("data")
    
    with pytest.raises(FileExistsError):
        make_blast_db(str(faa1), "db_exists", overwrite=False)
    
    # Reset mock stdout for next call
    mock_popen_instance.stdout = io.StringIO("building...")
    # This should now hit the 'except FileNotFoundError' branch for other extensions
    make_blast_db(str(faa1), "db_exists", overwrite=True)
    assert os.path.exists("db_exists.psq") == False

    # 5. Missing binary
    mock_run.side_effect = FileNotFoundError
    with pytest.raises(ValueError):
        make_blast_db(str(faa1), "db")
    mock_run.side_effect = None

    # 6. Bad extension
    with pytest.raises(ValueError):
        make_blast_db(str(bad_ext), "db")

    # 7. Check if temporary file is removed
    # (Checking for topiary_makeblastdb_*.faa)
    import glob
    assert len(glob.glob("topiary_makeblastdb_*.faa")) == 0

