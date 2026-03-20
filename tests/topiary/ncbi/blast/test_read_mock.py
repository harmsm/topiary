
import pytest
import os
from topiary.ncbi.blast.read import read_blast_xml, check_for_cpu_limit

def test_read_blast_xml_edge_cases(tmpdir):
    
    # 1. Directory with no .xml (line 217)
    empty_dir = tmpdir.mkdir("empty")
    with pytest.raises(ValueError):
        read_blast_xml(str(empty_dir))
        
    # 2. Invalid xml_input type (line 231-233)
    with pytest.raises(ValueError):
        read_blast_xml(123)
        
    # 3. List with missing file (line 239-240)
    with pytest.raises(ValueError):
        read_blast_xml(["missing.xml"])

def test_check_for_cpu_limit_not_found():
    with pytest.raises(FileNotFoundError):
        check_for_cpu_limit("missing.xml")
