import pytest

from topiary import io
import numpy as np
import pandas as pd

import warnings, os, shutil, re

def test_read_dataframe(dataframe_good_files,test_dataframes):
    """
    Test read dataframe function.
    """

    ref_df = test_dataframes["good-df"]

    for f in dataframe_good_files:

        # Print f in case parsing dies... so we know which file causes failure.
        print(f)

        # Read file and make sure it does not throw warning.
        with warnings.catch_warnings():
            warnings.simplefilter("error")

            df = io.read_dataframe(f)
            assert len(df) == len(ref_df)

            # Make sure expected columns are present
            df.uid
            df.keep

    # Check reading a dataframe
    df = io.read_dataframe(ref_df)
    assert len(df) == len(ref_df)
    assert df is not ref_df # make sure it's a copy
    df.uid
    df.keep

    # Make sure dies with useful error
    bad_inputs = [1,-1,1.5,None,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            io.read_dataframe(b)

    # Make sure raises file not found if a file is not passed
    with pytest.raises(FileNotFoundError):
        io.read_dataframe("not_really_a_file.txt")

def test_write_dataframe(test_dataframes,tmpdir):

    df = test_dataframes["good-df"]

    bad_df = [pd.DataFrame,pd.DataFrame({"test":[1]}),None,1,"string",str]
    for b in bad_df:
        with pytest.raises(ValueError):
            io.write_dataframe(b,"output_file.csv")

    bad_out_file = [pd.DataFrame,pd.DataFrame({"test":[1]}),None,1,str]
    for b in bad_out_file:
        with pytest.raises(ValueError):
            io.write_dataframe(df,b)

    def _check_written_out(df,out,sep):
        assert os.path.isfile(out)
        with open(out) as f:
            for line in f:
                assert len(line.split(sep)) == len(df.columns)

    out = os.path.join(tmpdir,"stupid.csv")
    f = open(out,"w")
    f.write("stuff")
    f.close()
    with pytest.raises(FileExistsError):
        io.write_dataframe(df,out_file=out)

    io.write_dataframe(df,out_file=out,overwrite=True)
    _check_written_out(df,out,",")

    # Write out as csv with non-standard extension
    out = os.path.join(tmpdir,"some_file.txt")
    io.write_dataframe(df,out_file=out)
    _check_written_out(df,out,",")

    out = os.path.join(tmpdir,"some_file.tsv")
    io.write_dataframe(df,out_file=out)
    _check_written_out(df,out,"\t")

    out = os.path.join(tmpdir,"some_file.xlsx")
    io.write_dataframe(df,out_file=out)
    assert os.path.exists(out)

def _validate_seq_writer():
    """
    This function is tested implicitly by test_write_fasta and test_write_phy.
    I'm using those unit tests because I can validate written output.
    """

    pass

def test_write_fasta(test_dataframes,tmpdir):

    def _check_output_file(out_file,num_columns,check_length=None):

        f = open(out_file)
        lines = f.readlines()
        f.close()

        has_gap = False
        for l in lines:
            if l.startswith(">"):
                col = l[1:].split("|")
                assert len(col) == num_columns
                assert len(col[0]) == 10
            else:
                if re.search("-",l):
                    has_gap = True

        if check_length is not None:
            assert len(lines) == check_length*2

        return has_gap

    df = test_dataframes["good-df_with-good-alignment"]

    # Basic read/write
    out_file = os.path.join(tmpdir,"output.fasta")
    io.write_fasta(df,out_file)
    _check_output_file(out_file,num_columns=3)

    # -------------------------------------------------------------------------
    # output files

    bad_output_files = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2)]
    for b in bad_output_files:
        with pytest.raises(ValueError):
            io.write_fasta(df,b)

    # -------------------------------------------------------------------------
    # seq_column

    os.remove(out_file)
    io.write_fasta(df,out_file=out_file)
    has_gaps = _check_output_file(out_file,num_columns=3)
    assert not has_gaps

    os.remove(out_file)
    io.write_fasta(df,out_file=out_file,seq_column="alignment")
    has_gaps = _check_output_file(out_file,num_columns=3)
    assert has_gaps

    bad_seq_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2),"not_in_df"]
    for b in bad_seq_columns:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,seq_column=b)

    # -------------------------------------------------------------------------
    # columns

    bad_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2)]
    for b in bad_columns:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,label_columns=b)

    os.remove(out_file)
    io.write_fasta(df,out_file,label_columns=["start"])
    has_gaps = _check_output_file(out_file,num_columns=2)

    # -------------------------------------------------------------------------
    # write_only_keepers

    bad_write_only = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_write_only:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,write_only_keepers=b)

    no_keep_df = df.copy()
    no_keep_df.loc[1:,"keep"] = False
    os.remove(out_file)
    io.write_fasta(no_keep_df,out_file,write_only_keepers=True)
    has_gaps = _check_output_file(out_file,num_columns=3,check_length=1)

    no_keep_df.keep = False
    os.remove(out_file)
    with pytest.raises(RuntimeError):
        io.write_fasta(no_keep_df,out_file,write_only_keepers=True)

    # -------------------------------------------------------------------------
    # empty_char

    bad_empty = [1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                 str,(1,2),True]
    for b in bad_empty:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,empty_char=b)

    # Should throw runtime error because sequences all look empty!
    with pytest.raises(RuntimeError):
        io.write_fasta(df,out_file=out_file,empty_char="ACDEFGHIKLMNPQRSTVWYZ-")

    # -------------------------------------------------------------------------
    # clean_sequence

    bad_clean = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_clean:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,clean_sequence=b)

    to_clean_df = df.copy()
    to_clean_df.loc[:,"sequence"] = "STUPIDX?STUPID"
    to_clean_df.loc[:,"length"] = len("STUPIDX?STUPID")
    io.write_fasta(to_clean_df,out_file=out_file,clean_sequence=True)
    f = open(out_file)
    lines = f.readlines()
    f.close()
    assert lines[1].strip() == "ST-PID--ST-PID"

    # -------------------------------------------------------------------------
    # overwrite

    bad_overwrite = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_overwrite:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,overwrite=b)

    os.remove(out_file)
    io.write_fasta(df,out_file)
    # Should throw error because overwrite is False by default
    with pytest.raises(FileExistsError):
        io.write_fasta(df,out_file)

    # Should throw error
    with pytest.raises(FileExistsError):
        io.write_fasta(df,out_file,overwrite=False)

    # Should work because we have the file in place
    io.write_fasta(df,out_file,overwrite=True)

def test_write_phy(test_dataframes,tmpdir):

    def _check_output_file(out_file,check_length=None):

        f = open(out_file)
        lines = f.readlines()
        f.close()

        header = lines[0].split()
        num_seqs = int(header[0])
        seq_length = int(header[1])

        if check_length is not None:
            assert num_seqs == check_length

        counter = 0
        has_gap = False
        for l in lines[2:]:

            if counter == 0:
                assert len(l.strip()) == 10
                counter += 1
            else:
                assert len(l.strip()) == seq_length
                if re.search("-",l):
                    has_gap = True

                counter = 0

        return has_gap

    df = test_dataframes["good-df_with-good-alignment"]

    # Basic read/write
    out_file = os.path.join(tmpdir,"output.fasta")
    io.write_phy(df,out_file)
    _check_output_file(out_file)

    # -------------------------------------------------------------------------
    # output files

    bad_output_files = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2)]
    for b in bad_output_files:
        with pytest.raises(ValueError):
            io.write_phy(df,b)

    # -------------------------------------------------------------------------
    # seq_column

    os.remove(out_file)
    io.write_phy(df,out_file=out_file,seq_column="sequence")
    has_gaps = _check_output_file(out_file)
    assert not has_gaps

    os.remove(out_file)
    io.write_phy(df,out_file=out_file,seq_column="alignment")
    has_gaps = _check_output_file(out_file)
    assert has_gaps

    bad_seq_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2),"not_in_df"]
    for b in bad_seq_columns:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,seq_column=b)

    # Send in alignment column with a short sequence (not all same length)
    bad_align_df = df.copy()
    bad_align_df.loc[0,"alignment"] = "MTG"
    os.remove(out_file)
    with pytest.raises(ValueError):
        io.write_phy(bad_align_df,out_file=out_file,seq_column="alignment")

    # -------------------------------------------------------------------------
    # write_only_keepers

    bad_write_only = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_write_only:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,write_only_keepers=b)

    no_keep_df = df.copy()
    no_keep_df.loc[1:,"keep"] = False
    io.write_phy(no_keep_df,out_file,write_only_keepers=True)
    has_gaps = _check_output_file(out_file,check_length=1)

    no_keep_df.keep = False
    os.remove(out_file)
    with pytest.raises(ValueError):
        io.write_phy(no_keep_df,out_file,write_only_keepers=True)

    # -------------------------------------------------------------------------
    # empty_char

    bad_empty = [1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                 str,(1,2),True]
    for b in bad_empty:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,empty_char=b)

    # Should throw runtime error because sequences all look empty!
    with pytest.raises(RuntimeError):
        io.write_phy(df,out_file=out_file,empty_char="ACDEFGHIKLMNPQRSTVWYZ-")

    # -------------------------------------------------------------------------
    # clean_sequence

    bad_clean = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_clean:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,clean_sequence=b)

    to_clean_df = df.copy()
    to_clean_df.loc[0,"alignment"] = "?LUFLFF---"
    io.write_phy(to_clean_df,out_file=out_file,clean_sequence=True)
    f = open(out_file)
    lines = f.readlines()
    f.close()
    assert lines[3].strip() == "-L-FLFF---"

    # -------------------------------------------------------------------------
    # overwrite

    bad_overwrite = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_overwrite:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,overwrite=b)

    os.remove(out_file)
    io.write_phy(df,out_file)
    # Should throw error because overwrite is False by default
    with pytest.raises(FileExistsError):
        io.write_phy(df,out_file)

    # Should throw error
    with pytest.raises(FileExistsError):
        io.write_phy(df,out_file,overwrite=False)

    # Should work because we have the file in place
    io.write_phy(df,out_file,overwrite=True)

def test_ncbi_blast_xml_to_df(xml,tmpdir):

    # Pass in a single xml file, not in a list
    xml_file = xml["good.xml"]
    df = io.ncbi_blast_xml_to_df(xml_file)
    assert len(df) == 19

    # Pass two xml files (indetical, so should end up with single 19-length df)
    df = io.ncbi_blast_xml_to_df([xml_file,xml_file])
    assert len(df) == 19

    # Pass directory with an xml file
    xml_files_dir = os.path.join(tmpdir,"xml_files")
    os.mkdir(xml_files_dir)
    shutil.copy(xml_file,
                os.path.join(xml_files_dir,os.path.split(xml_file)[-1]))
    df = io.ncbi_blast_xml_to_df(xml_files_dir)
    assert len(df) == 19

    # Passing a stupid xml file with broken format ... not a great test, but
    # should at least throw error of some sort.
    xml_file = xml["bad.xml"]
    with pytest.raises(ValueError):
        df = io.ncbi_blast_xml_to_df(xml_file)
