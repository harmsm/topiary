import pytest

import sys

import topiary
from topiary.cli_scripts.alignment_to_ancestors import main
from topiary.pipeline import alignment_to_ancestors


def test_main_wraps_the_pipeline_function(mocker):

    wrap = mocker.patch("topiary.cli_scripts.alignment_to_ancestors.wrap_function")

    main(argv=["--df","dataframe.csv"])

    wrap.assert_called_once()
    assert wrap.call_args[0][0] is alignment_to_ancestors
    assert wrap.call_args.kwargs["argv"] == ["--df","dataframe.csv"]
    assert wrap.call_args.kwargs["optional_arg_types"] == {}


def test_main_defaults_to_sys_argv(mocker):

    wrap = mocker.patch("topiary.cli_scripts.alignment_to_ancestors.wrap_function")
    mocker.patch.object(sys,"argv",["topiary-alignment-to-ancestors","--df","x.csv"])

    main()

    assert wrap.call_args.kwargs["argv"] == ["--df","x.csv"]
