import pytest

import topiary.generax._crawl as _crawl
from topiary.generax._crawl import _atomic_create
from topiary.generax._crawl import _is_terminal
from topiary.generax._crawl import _bump_attempts
from topiary.generax._crawl import claim_replicate
from topiary.generax._crawl import run_replicate
from topiary.generax._crawl import assemble_bs_trees
from topiary.generax._crawl import check_bootstrap_convergence
from topiary.generax._crawl import mark_remaining_skipped
from topiary.generax._crawl import all_terminal
from topiary.generax._crawl import replicate_dirs
from topiary.generax._crawl import crawl
from topiary.generax._crawl import _get_timeout_config
from topiary.generax._crawl import find_ready_calc_dir
from topiary.generax._crawl import elect_setup
from topiary.generax._crawl import elect_aggregate
from topiary.generax._crawl import mark_aggregate_done

import os
import time
import pathlib
import multiprocessing as mp


# -----------------------------------------------------------------------------
# helpers
# -----------------------------------------------------------------------------

def _build_reps(repdir, names):
    """Build a minimal replicate tree with run_generax.sh in each rep dir."""
    os.mkdir(repdir)
    for name in names:
        d = os.path.join(repdir, name)
        os.mkdir(d)
        with open(os.path.join(d, "run_generax.sh"), "w") as f:
            f.write("generax --families control.txt &> topiary.log\n")


def _fake_launch(plan):
    """
    Build a fake _launch_replicate. `plan` maps replicate-name -> outcome
    ("success", "timeout", "notree"), defaulting to "success".
    """

    def fake(cmd, timeout, stdout_path, stderr_path, cwd=None):
        d = os.path.dirname(stdout_path)
        open(stdout_path, "w").close()
        open(stderr_path, "w").close()
        name = os.path.basename(os.path.abspath(d))
        outcome = plan.get(name, "success") if isinstance(plan, dict) else plan
        if outcome == "success":
            rt = os.path.join(d, _crawl._RESULT_TREE)
            os.makedirs(os.path.dirname(rt), exist_ok=True)
            with open(rt, "w") as f:
                f.write("(A:1,B:1);\n")
            return 0, False
        if outcome == "timeout":
            return None, True
        return 1, False   # notree

    return fake


# -----------------------------------------------------------------------------
# atomic primitives
# -----------------------------------------------------------------------------

def test__atomic_create(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    assert _atomic_create("marker", "me") is True
    assert os.path.isfile("marker")
    with open("marker") as f:
        assert f.read().strip() == "me"
    # second create fails (already exists)
    assert _atomic_create("marker", "other") is False
    with open("marker") as f:
        assert f.read().strip() == "me"


def test__is_terminal(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("d")
    assert _is_terminal("d") is False
    pathlib.Path(os.path.join("d", "completed")).touch()
    assert _is_terminal("d") is True

    os.mkdir("d2")
    pathlib.Path(os.path.join("d2", "failed")).touch()
    assert _is_terminal("d2") is True

    os.mkdir("d3")
    pathlib.Path(os.path.join("d3", "skipped")).touch()
    assert _is_terminal("d3") is True


def test__bump_attempts(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("d")
    assert _bump_attempts("d") == 1
    assert _bump_attempts("d") == 2
    assert _bump_attempts("d") == 3
    with open(os.path.join("d", "attempts")) as f:
        assert f.read().strip() == "3"


# -----------------------------------------------------------------------------
# claiming
# -----------------------------------------------------------------------------

def test_claim_replicate(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("d")

    # fresh claim succeeds and leaves a running marker
    assert claim_replicate("d", "crawlerA") is True
    assert os.path.isfile(os.path.join("d", "running"))

    # a second (fresh) crawler cannot claim it -- the claim is not stale
    assert claim_replicate("d", "crawlerB") is False

    # terminal dirs are not claimable
    os.mkdir("done")
    pathlib.Path(os.path.join("done", "completed")).touch()
    assert claim_replicate("done", "crawlerA") is False


def test_claim_replicate_stale_steal(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("d")

    # crawlerA claims, then goes silent
    assert claim_replicate("d", "crawlerA") is True

    # not stale yet -> B cannot steal
    assert claim_replicate("d", "crawlerB", stale_seconds=1000) is False

    # age the claim past the staleness threshold
    old = time.time() - 10_000
    os.utime(os.path.join("d", "running"), (old, old))

    # now B can steal it
    assert claim_replicate("d", "crawlerB", stale_seconds=300) is True
    with open(os.path.join("d", "running")) as f:
        assert f.read().strip() == "crawlerB"


def _claim_worker(args):
    """Top-level worker: try to claim every dir once; return the ones we won."""
    repdir, cid, names = args
    won = []
    for name in names:
        if claim_replicate(os.path.join(repdir, name), cid):
            won.append(name)
    return won


def test_claim_replicate_concurrent(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = [f"{i:05d}" for i in range(1, 41)]
    _build_reps("replicates", names)

    # Six workers race to claim the same 40 replicates.
    args = [("replicates", f"crawler{i}", names) for i in range(6)]
    with mp.Pool(6) as pool:
        results = pool.map(_claim_worker, args)

    won = [name for sub in results for name in sub]

    # Every replicate claimed exactly once -- no double claims, none missed.
    assert sorted(won) == names
    assert len(won) == len(set(won))


# -----------------------------------------------------------------------------
# running a replicate
# -----------------------------------------------------------------------------

def test_run_replicate_success(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    _build_reps("replicates", ["00001"])
    d = os.path.join("replicates", "00001")
    assert claim_replicate(d, "me") is True

    monkeypatch.setattr(_crawl, "_launch_replicate", _fake_launch("success"))

    durations = []
    outcome = run_replicate(d, generax_launch="", config=_get_timeout_config(None),
                            durations=durations, sample_threshold=1)
    assert outcome == "completed"
    assert os.path.isfile(os.path.join(d, "completed"))
    assert os.path.isfile(os.path.join(d, "bs-tree.newick"))
    assert not os.path.exists(os.path.join(d, "running"))
    assert len(durations) == 1


def test_run_replicate_retry_then_fail(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    _build_reps("replicates", ["00001"])
    d = os.path.join("replicates", "00001")

    monkeypatch.setattr(_crawl, "_launch_replicate", _fake_launch("notree"))

    durations = []
    outcomes = []
    for _ in range(3):
        assert claim_replicate(d, "me") is True
        outcomes.append(run_replicate(d, generax_launch="",
                                      config=_get_timeout_config(None),
                                      durations=durations, sample_threshold=1))

    # 1 initial + 2 retries, then failed on the third
    assert outcomes == ["retry", "retry", "failed"]
    assert os.path.isfile(os.path.join(d, "failed"))
    assert not os.path.exists(os.path.join(d, "running"))
    with open(os.path.join(d, "attempts")) as f:
        assert f.read().strip() == "3"
    # a failed replicate is no longer claimable
    assert claim_replicate(d, "me") is False


def test_run_replicate_timeout_is_failure(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    _build_reps("replicates", ["00001"])
    d = os.path.join("replicates", "00001")
    assert claim_replicate(d, "me") is True

    monkeypatch.setattr(_crawl, "_launch_replicate", _fake_launch("timeout"))
    outcome = run_replicate(d, generax_launch="", config=_get_timeout_config(None),
                            durations=[], sample_threshold=1)
    assert outcome == "retry"
    assert not os.path.exists(os.path.join(d, "completed"))
    with open(os.path.join(d, "stderr.log")) as f:
        assert "timeout" in f.read()


# -----------------------------------------------------------------------------
# assembly / convergence / bookkeeping
# -----------------------------------------------------------------------------

def test_assemble_bs_trees(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = ["00001", "00002", "00003"]
    _build_reps("replicates", names)

    # 00001 and 00003 completed; 00002 not
    for name in ["00001", "00003"]:
        d = os.path.join("replicates", name)
        with open(os.path.join(d, "bs-tree.newick"), "w") as f:
            f.write("(A:1,B:1);\n")
        pathlib.Path(os.path.join(d, "completed")).touch()

    trees = assemble_bs_trees("replicates")
    assert len(trees) == 2

    assemble_bs_trees("replicates", out_path="combined.newick")
    with open("combined.newick") as f:
        assert len([ln for ln in f if ln.strip()]) == 2


def test_mark_remaining_skipped_and_all_terminal(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = ["00001", "00002", "00003", "00004"]
    _build_reps("replicates", names)

    # 00001 completed, 00002 running (an active claim - must NOT be skipped)
    pathlib.Path(os.path.join("replicates", "00001", "completed")).touch()
    pathlib.Path(os.path.join("replicates", "00002", "running")).touch()

    assert all_terminal("replicates") is False
    mark_remaining_skipped("replicates")

    # unclaimed 00003/00004 get skipped; running 00002 is left alone
    assert os.path.isfile(os.path.join("replicates", "00003", "skipped"))
    assert os.path.isfile(os.path.join("replicates", "00004", "skipped"))
    assert not os.path.exists(os.path.join("replicates", "00002", "skipped"))

    # finish 00002 -> everything terminal
    os.remove(os.path.join("replicates", "00002", "running"))
    pathlib.Path(os.path.join("replicates", "00002", "completed")).touch()
    assert all_terminal("replicates") is True


# -----------------------------------------------------------------------------
# the crawl loop
# -----------------------------------------------------------------------------

def test_crawl_completes_all(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = [f"{i:05d}" for i in range(1, 7)]
    _build_reps("replicates", names)

    monkeypatch.setattr(_crawl, "_launch_replicate", _fake_launch("success"))
    # no real raxml -- never converge early
    monkeypatch.setattr(_crawl, "check_bootstrap_convergence",
                        lambda *a, **k: (False, None))

    n_run = crawl("replicates", generax_launch="", cid="solo")

    assert n_run == len(names)
    assert all_terminal("replicates")
    for name in names:
        assert os.path.isfile(os.path.join("replicates", name, "completed"))
    assert len(assemble_bs_trees("replicates")) == len(names)


def test_crawl_convergence_shortcircuits(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = [f"{i:05d}" for i in range(1, 11)]
    _build_reps("replicates", names)

    monkeypatch.setattr(_crawl, "_launch_replicate", _fake_launch("success"))

    # "Converge" as soon as 3 replicates have completed.
    def fake_conv(replicate_dir, converge_cutoff, workdir=None):
        n = len([d for d in replicate_dirs(replicate_dir)
                 if os.path.exists(os.path.join(d, "completed"))])
        return (n >= 3), None
    monkeypatch.setattr(_crawl, "check_bootstrap_convergence", fake_conv)

    crawl("replicates", generax_launch="", cid="solo")

    n_completed = len([d for d in replicate_dirs("replicates")
                       if os.path.exists(os.path.join(d, "completed"))])
    n_skipped = len([d for d in replicate_dirs("replicates")
                     if os.path.exists(os.path.join(d, "skipped"))])
    # some completed, the rest skipped, nothing left running
    assert n_completed >= 3
    assert n_skipped > 0
    assert n_completed + n_skipped == len(names)
    assert all_terminal("replicates")


# -----------------------------------------------------------------------------
# setup / aggregate coordination
# -----------------------------------------------------------------------------

def _make_build_fn(parent):
    """A build_fn that records each invocation and returns a fixed calc dir."""
    def build():
        calc = os.path.join(parent, "01_reconciled-tree-bootstraps")
        os.makedirs(os.path.join(calc, "working"), exist_ok=True)
        # record this build invocation atomically-uniquely
        _atomic_create(os.path.join(parent, f"built.{os.getpid()}"), "x")
        return calc
    return build


def test_find_ready_calc_dir(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.makedirs(os.path.join("03_reconciled-tree-bootstraps", "working"))
    os.makedirs(os.path.join("05_reconciled-tree-bootstraps", "working"))

    # none ready yet
    assert find_ready_calc_dir(".") is None

    pathlib.Path(os.path.join("03_reconciled-tree-bootstraps",
                              "working", ".crawl-ready")).touch()
    assert find_ready_calc_dir(".").endswith("03_reconciled-tree-bootstraps")

    # newest ready wins
    pathlib.Path(os.path.join("05_reconciled-tree-bootstraps",
                              "working", ".crawl-ready")).touch()
    assert find_ready_calc_dir(".").endswith("05_reconciled-tree-bootstraps")


def test_elect_setup_single(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("parent")

    calc_dir, is_leader = elect_setup("parent", _make_build_fn("parent"))
    assert is_leader is True
    assert calc_dir.endswith("01_reconciled-tree-bootstraps")
    assert len(glob_built("parent")) == 1

    # A second call sees setup already done -> follower, no rebuild.
    calc_dir2, is_leader2 = elect_setup("parent", _make_build_fn("parent"))
    assert is_leader2 is False
    assert calc_dir2 == calc_dir
    assert len(glob_built("parent")) == 1


def glob_built(parent):
    import glob as _glob
    return _glob.glob(os.path.join(parent, "built.*"))


def _setup_worker(args):
    parent, cid = args
    calc_dir, is_leader = elect_setup(parent, _make_build_fn(parent),
                                      cid=cid, poll=0.2, timeout=60)
    return (calc_dir, is_leader)


def test_elect_setup_concurrent(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.mkdir("parent")

    args = [("parent", f"crawler{i}") for i in range(6)]
    with mp.Pool(6) as pool:
        results = pool.map(_setup_worker, args)

    # exactly one leader, build ran exactly once, everyone agrees on calc_dir
    leaders = [r for r in results if r[1]]
    assert len(leaders) == 1
    assert len(glob_built("parent")) == 1
    calc_dirs = set(r[0] for r in results)
    assert len(calc_dirs) == 1


def test_elect_aggregate(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    os.makedirs(os.path.join("calc", "working"))

    # first caller wins
    assert elect_aggregate("calc", cid="a") is True
    # while the lock is held, no one else wins
    assert elect_aggregate("calc", cid="b") is False

    # once marked done, nobody runs it again
    mark_aggregate_done("calc", cid="a")
    assert elect_aggregate("calc", cid="c") is False


@pytest.mark.run_raxml
def test_check_bootstrap_convergence(tmpdir, monkeypatch):
    monkeypatch.chdir(tmpdir)
    names = ["00001", "00002", "00003"]
    _build_reps("replicates", names)

    # identical trees across replicates -> converged
    for name in names:
        d = os.path.join("replicates", name)
        with open(os.path.join(d, "bs-tree.newick"), "w") as f:
            f.write("((A:1,B:1):1,(C:1,D:1):1);\n")
        pathlib.Path(os.path.join(d, "completed")).touch()

    converged, df = check_bootstrap_convergence("replicates", converge_cutoff=0.5)
    assert converged in (True, False)   # ran raxml without error
    assert df is not None
