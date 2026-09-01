import pytest

import sys

import topiary
from topiary.cli_scripts.create_report import main
from topiary.reports import tree_report


def test_main_wraps_tree_report(mocker):

    wrap = mocker.patch("topiary.cli_scripts.create_report.wrap_function")

    main(argv=["--output_directory","out"])

    wrap.assert_called_once()
    assert wrap.call_args[0][0] is tree_report
    assert wrap.call_args.kwargs["argv"] == ["--output_directory","out"]

    # ancestor_directory defaults to None in tree_report, so the CLI has to
    # declare its type explicitly or argparse cannot build the flag.
    assert wrap.call_args.kwargs["optional_arg_types"] == {"ancestor_directory":str}


def test_main_defaults_to_sys_argv(mocker):

    wrap = mocker.patch("topiary.cli_scripts.create_report.wrap_function")
    mocker.patch.object(sys,"argv",["topiary-create-report","--output_directory","o"])

    main()

    assert wrap.call_args.kwargs["argv"] == ["--output_directory","o"]
