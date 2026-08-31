import pytest

from topiary.pipeline.bootstrap_reconcile import bootstrap_reconcile
from topiary.pipeline.bootstrap_reconcile import _input_bootstrap_dir
from topiary._private.interface import WrappedFunctionException

import os


def _fake_flow(mocker, input_status="complete",
               is_leader=True, all_terminal=True, is_aggregator=True):
    """Mock out the crawler engine + software validation for the pipeline flow."""
    mocker.patch("topiary.pipeline.bootstrap_reconcile.installed.validate_stack")

    sv = mocker.Mock()
    sv.status = input_status
    mocker.patch("topiary.pipeline.bootstrap_reconcile.Supervisor", return_value=sv)

    mocker.patch("topiary.generax._crawl.elect_setup",
                 return_value=("06_reconciled-tree-bootstraps", is_leader))
    crawl = mocker.patch("topiary.generax._crawl.crawl")
    mocker.patch("topiary.generax._crawl.all_terminal", return_value=all_terminal)
    mocker.patch("topiary.generax._crawl.elect_aggregate", return_value=is_aggregator)
    aggregate = mocker.patch("topiary.generax._crawl.aggregate_bootstrap")
    mocker.patch("topiary.generax._crawl.mark_aggregate_done")
    report = mocker.patch("topiary.reports.pipeline_report")

    return {"crawl": crawl, "aggregate": aggregate, "report": report}


def _make_run_dir(tmpdir, input_name="05_gene-tree-bootstraps"):
    pdir = os.path.join(tmpdir, "run")
    os.makedirs(os.path.join(pdir, input_name))
    return pdir


def test_bootstrap_reconcile_flow(tmpdir, mocker):

    mocks = _fake_flow(mocker)
    pdir = _make_run_dir(tmpdir)

    bootstrap_reconcile(pdir)

    # crawler ran, aggregation ran, report written
    mocks["crawl"].assert_called_once()
    mocks["aggregate"].assert_called_once()
    mocks["report"].assert_called_once()

    # the launch prefix is threaded into crawl
    assert mocks["crawl"].call_args.kwargs["generax_launch"] == ""


def test_bootstrap_reconcile_generax_launch(tmpdir, mocker):

    mocks = _fake_flow(mocker)
    pdir = _make_run_dir(tmpdir)

    bootstrap_reconcile(pdir, generax_launch="mpirun -np 8")
    assert mocks["crawl"].call_args.kwargs["generax_launch"] == "mpirun -np 8"


def test_bootstrap_reconcile_follower_no_aggregate(tmpdir, mocker):

    # A crawler that is neither leader nor the elected aggregator crawls but does
    # not aggregate or write the report.
    mocks = _fake_flow(mocker, is_leader=False, is_aggregator=False)
    pdir = _make_run_dir(tmpdir)

    bootstrap_reconcile(pdir)
    mocks["crawl"].assert_called_once()
    mocks["aggregate"].assert_not_called()
    mocks["report"].assert_not_called()


def test_bootstrap_reconcile_no_aggregate_until_terminal(tmpdir, mocker):

    # If replicates are not all terminal, this crawler does not aggregate.
    mocks = _fake_flow(mocker, all_terminal=False)
    pdir = _make_run_dir(tmpdir)

    bootstrap_reconcile(pdir)
    mocks["crawl"].assert_called_once()
    mocks["aggregate"].assert_not_called()


def test_bootstrap_reconcile_errors(tmpdir, mocker):

    # previous_run_dir does not exist
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile(os.path.join(tmpdir, "nope"))

    # previous_run_dir exists but has no input bootstraps directory
    mocker.patch("topiary.pipeline.bootstrap_reconcile.installed.validate_stack")
    empty = os.path.join(tmpdir, "empty")
    os.mkdir(empty)
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile(empty)

    # input calculation not complete
    _fake_flow(mocker, input_status="running")
    pdir = _make_run_dir(tmpdir)
    with pytest.raises(WrappedFunctionException):
        bootstrap_reconcile(pdir)


def test__input_bootstrap_dir(tmpdir, monkeypatch):

    monkeypatch.chdir(tmpdir)
    os.mkdir("03_gene-tree-bootstraps")
    os.mkdir("05_gene-tree-bootstraps")
    os.mkdir("06_reconciled-tree-bootstraps")   # must be ignored as input

    num, name = _input_bootstrap_dir(".")
    assert num == 5
    assert name == "05_gene-tree-bootstraps"
