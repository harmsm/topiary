#!/usr/bin/env python3
"""
Check for installed packages on command line.
"""

import topiary

import os

def main(argv=None):

    required = topiary._private.software_requirements

    print(70*"=")
    print("Checking for installed external software")
    print(70*"=")
    print("",flush=True)

    to_check = []
    for r in required:
        to_check.append({"program":r,
                         "min_version":required[r],
                         "must_pass":False})

    try:
        topiary._private.installed.validate_stack(to_check)
    except RuntimeError as e:

        # validate_stack raises a RuntimeError whose message already
        # distinguishes between programs that are missing from the $PATH,
        # programs whose version is too low, and programs that are present but
        # crash when run (e.g. a binary compiled for a different architecture).
        # Surface that message directly rather than assuming a $PATH problem.
        print(str(e))

        print("\nThe current $PATH visible to topiary contains the following")
        print("directories:")

        path_dirs = os.environ["PATH"].split(os.pathsep)
        for d in path_dirs:
            print(f"    {d}")
        print("")

if __name__ == "__main__":
    main()
