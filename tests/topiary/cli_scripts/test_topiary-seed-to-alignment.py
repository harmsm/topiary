import pytest

import sys

import topiary
from topiary.cli_scripts.seed_to_alignment import main
from topiary.pipeline import seed_to_alignment


def test_main_wraps_the_pipeline_function(mocker):

    wrap = mocker.patch("topiary.cli_scripts.seed_to_alignment.wrap_function")

    main(argv=["--seed_df","seed.csv"])

    wrap.assert_called_once()

    # The console script must hand wrap_function the pipeline entry point, not
    # something else -- this is the only thing connecting `topiary-seed-to-
    # alignment` on the command line to the code that does the work.
    assert wrap.call_args[0][0] is seed_to_alignment
    assert wrap.call_args.kwargs["argv"] == ["--seed_df","seed.csv"]
    assert wrap.call_args.kwargs["optional_arg_types"] == {}


def test_main_defaults_to_sys_argv(mocker):

    wrap = mocker.patch("topiary.cli_scripts.seed_to_alignment.wrap_function")
    mocker.patch.object(sys,"argv",["topiary-seed-to-alignment","--seed_df","x.csv"])

    main()

    assert wrap.call_args.kwargs["argv"] == ["--seed_df","x.csv"]
