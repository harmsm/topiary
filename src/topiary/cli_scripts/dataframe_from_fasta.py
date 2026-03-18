#!/usr/bin/env python3
"""
Create a topiary dataframe from a fasta file.
"""

import topiary
from topiary._private import check
from topiary._private.wrap import wrap_function

import pandas as pd
from Bio import SeqIO

import sys

def dataframe_from_fasta(fasta,
                         out="topiary-dataframe.csv",
                         description="unknown",
                         species="unknown",
                         overwrite=False):
    """
    Create a topiary dataframe from a fasta file.

    Parameters
    ----------
    fasta : str
        fasta file to read.
    out : str, default="topiary-dataframe.csv"
        output file for the dataframe.
    description : str, default="unknown"
        placeholder for 'name' column in the dataframe if not parsed from fasta.
    species : str, default="unknown"
        placeholder for 'species' column in the dataframe.
    overwrite : bool, default=False
        whether or not to overwrite an existing output file.

    Returns
    -------
    pandas.DataFrame
        topiary dataframe.
    """

    # Create lists to hold data
    names = []
    sequences = []

    # Read fasta file
    for record in SeqIO.parse(fasta, "fasta"):
        names.append(record.description)
        sequences.append(str(record.seq))

    # Create dataframe
    df_dict = {
        "name": names,
        "species": [species for _ in range(len(names))],
        "sequence": sequences,
        "alignment": sequences,
        "keep": [True for _ in range(len(names))]
    }

    df = pd.DataFrame(df_dict)

    # Validate dataframe (will add uids etc.)
    df = check.check_topiary_dataframe(df)

    return df

def main(argv=None):
    """
    Main function for the script.
    """

    if argv is None:
        argv = sys.argv[1:]

    # Set arg types for args with None as default
    optional_arg_types = {"overwrite": bool}

    extra_args = [("out", {"type": str})]

    description = \
    """
    Create a topiary dataframe from a fasta file. This script generates a
    minimal dataframe (uid, name, species, sequence, alignment, keep) that can
    be used to start a topiary pipeline. It does not perform any sophisticated
    parsing of the fasta headers.
    """

    # Wrap and run function
    ret, args = wrap_function(dataframe_from_fasta,
                              argv=argv,
                              optional_arg_types=optional_arg_types,
                              extra_args=extra_args,
                              description=description)

    out_file = args.__dict__["out"]
    overwrite = args.__dict__["overwrite"]
    topiary.write_dataframe(ret, out_file, overwrite=overwrite)

if __name__ == "__main__":
    main()
