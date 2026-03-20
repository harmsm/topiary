import pytest
import topiary
import numpy as np
import pandas as pd
import os, subprocess, sys

def test_main(tmpdir):

    # Create fasta file
    fasta_file = os.path.abspath(os.path.join(tmpdir, "test.fasta"))
    f = open(fasta_file, 'w')
    f.write(">seq1 name 1\nMAEGE\n")
    f.write(">seq2 name 2\nMAE-E\n")
    f.close()

    out_file = os.path.abspath(os.path.join(tmpdir, "output-df.csv"))

    # Get location of binary relative to python interpreter
    bin_dir = os.path.dirname(sys.executable)
    test_bin = os.path.join(bin_dir, "topiary-dataframe-from-fasta")

    if os.name == "nt":
        base_cmd = [sys.executable, test_bin]
    else:
        base_cmd = [test_bin]

    # Should run but fail because no arguments
    ret = subprocess.run(base_cmd)
    assert ret.returncode != 0

    # Test basic conversion
    cmd = base_cmd[:]
    cmd.extend([fasta_file, out_file])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.path.isfile(out_file)
    
    out_df = topiary.read_dataframe(out_file)
    assert len(out_df) == 2
    assert "uid" in out_df.columns
    assert np.array_equal(out_df.name, ["seq1 name 1", "seq2 name 2"])
    assert np.array_equal(out_df.sequence, ["MAEGE", "MAE-E"])
    assert np.array_equal(out_df.alignment, ["MAEGE", "MAE-E"])
    assert np.all(out_df.keep)
    assert np.all(out_df.species == "unknown")

    # Test with custom species and description
    os.remove(out_file)
    cmd = base_cmd[:]
    cmd.extend([fasta_file, out_file, "--species", "Homo sapiens", "--description", "my protein"])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    
    out_df = topiary.read_dataframe(out_file)
    assert np.all(out_df.species == "Homo sapiens")

    # Test overwrite fail (default is False)
    ret = subprocess.run(cmd)
    assert ret.returncode != 0

    # Test overwrite success
    cmd.append("--overwrite")
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
