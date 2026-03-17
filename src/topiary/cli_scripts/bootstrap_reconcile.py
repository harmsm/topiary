#!/usr/bin/env python3
"""
Command line interface to pipeline.bootstrap_reconcile
"""

import topiary
from topiary.pipeline import bootstrap_reconcile
from topiary._private.wrap import wrap_function

import sys

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Set arg types for args with None as default
    optional_arg_types = {"num_threads": int,
                          "threads_per_replicate": int}

    # Wrap and run function
    wrap_function(bootstrap_reconcile,
                  argv=argv,
                  description=bootstrap_reconcile.__doc__,
                  optional_arg_types=optional_arg_types)

if __name__ == "__main__":
    main()
