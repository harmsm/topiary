
import pytest
import topiary
import copy
import pandas as pd
import numpy as np

from topiary.ncbi._parse_ncbi_line import _grab_line_meta_data
from topiary.ncbi._parse_ncbi_line import parse_ncbi_line

def test__grab_line_meta_data(ncbi_lines):

    # This function is implicitly and rigorously tested by test_parse_ncbi_line.
    # I'd have to make separate test data manually, so I'll rely on
    # test_parse_ncbi_line. 
    # Just a sanity call to make sure it runs on a simple string.
    out = _grab_line_meta_data("crystal structure of a protein")
    assert out["structure"] is True
    assert out["isoform"] is False

def test_parse_ncbi_line(ncbi_lines):

    input_lines = ncbi_lines[0]
    ncbi_lines_parsed = ncbi_lines[1]

    for i, line in enumerate(input_lines):
        line_dict = parse_ncbi_line(line)

        for k in ncbi_lines_parsed[i]:
            assert line_dict[k] == ncbi_lines_parsed[i][k]

    # Test ability to grab other entry off the line rather than first entry
    line_dict = parse_ncbi_line(input_lines[0],accession="AIC51010")

    expected = copy.deepcopy(ncbi_lines_parsed[0])
    expected["line"] = "gb|AIC51010.1| LY96, partial [synthetic construct]"
    expected["accession"] = "AIC51010"
    expected["name"] = "LY96, partial"
    expected["species"] = "synthetic construct"
    expected["precursor"] = False
    expected["partial"] = True

    for k in expected:
        assert line_dict[k] == expected[k]

    # Edge cases for coverage
    
    # 1. pd.isnull(line)
    assert parse_ncbi_line(np.nan) is None
    assert parse_ncbi_line(None) is None
    assert parse_ncbi_line("") is None

    # 2. line.startswith(">")
    line_dict = parse_ncbi_line(">ref|NP_001186083.1| protein [Homo sapiens]")
    assert line_dict["accession"] == "NP_001186083"

    # 3. PDB format
    line_dict = parse_ncbi_line("pdb|1STN|A protein [Homo sapiens]")
    assert line_dict["accession"] == "1STN_A"

    # 4. Malformed line (IndexError)
    assert parse_ncbi_line("malformed line") is None

    # 5. Accession not in line (KeyError)
    line_dict = parse_ncbi_line("ref|NP_001186083.1| protein [Homo sapiens]", accession="mismatch")
    assert line_dict["accession"] == "mismatch"
    assert line_dict["line"] == "ref|NP_001186083.1| protein [Homo sapiens]"
