#!/usr/bin/env python
"""
Fail loudly if a coverage run collected nothing.

When coverage.py measures no files it does not error -- it emits a
CoverageWarning and then reports 0% for everything, so the only symptom is a
`--fail-under` failure claiming coverage collapsed. That reads like a huge
regression when the actual problem is that coverage was pointed at the wrong
copy of the code.

The way that happens here: `source = ["src/topiary"]` measures files *at that
path*, which is only where topiary lives when it is installed editable
(`pip install -e .`). CI installs non-editable, so the code that executes is in
site-packages, nothing under src/topiary is traced, and the total is 0%.
`source_pkgs = ["topiary"]` in pyproject.toml fixes that by following the
package wherever it is imported from -- this script exists to make a regression
of that kind say so.

Note that counting *measured files* is not enough to detect this. With a
directory-based `source`, coverage writes every file it finds under that
directory into the data file even when none of them were executed, so the file
count looks healthy while every file has zero executed lines. This checks
executed lines instead.

Usage
-----
    ./tests/check_coverage_collected.py [--min-files N]
"""

import argparse
import sys

import coverage


def main(argv=None):

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-files", type=int, default=50,
                        help="fail if fewer than this many files were measured "
                             "(default: 50)")

    args = parser.parse_args(argv)

    cov = coverage.Coverage()
    cov.load()
    data = cov.get_data()

    # A file only counts if lines in it actually ran. With a directory-based
    # `source`, coverage lists files it never traced, so len(measured_files())
    # stays high even when nothing was executed.
    with_lines = sorted(f for f in data.measured_files() if data.lines(f))
    executed = sum(len(data.lines(f)) for f in with_lines)

    if len(with_lines) >= args.min_files:
        print(f"coverage measured {len(with_lines)} files "
              f"({executed} executed lines)")
        return 0

    print(f"FAIL: coverage recorded executed lines in only "
          f"{len(with_lines)} files (expected at least {args.min_files}).")
    print(f"      {len(data.measured_files())} files are listed in the data "
          f"file, but {len(data.measured_files()) - len(with_lines)} of them "
          f"have no executed lines at all.")
    print()
    print("This almost always means coverage was watching a different copy of")
    print("topiary than the one the tests imported, not that coverage really")
    print("collapsed. Check:")
    print()
    print("  - Is topiary installed non-editable, so it runs from")
    print("    site-packages while coverage is pointed at src/topiary?")
    print("    `python -c 'import topiary; print(topiary.__file__)'` will say.")
    print("  - Does pyproject.toml still use source_pkgs = [\"topiary\"]")
    print("    rather than source = [\"src/topiary\"]?")
    print()

    if len(with_lines) > 0:
        print("Files that did record executed lines:")
        for path in with_lines[:20]:
            print(f"    {path}")

    return 1


if __name__ == "__main__":
    sys.exit(main())
