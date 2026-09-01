"""
Filesystem-coordinated crawler for gene/species tree reconciliation bootstrap
replicates.

Each bootstrap replicate is an independent generax run in its own directory.
Rather than orchestrating them with MPI/mpirun, topiary launches N independent
crawler processes (e.g. one per compute node via a SLURM job array). The
crawlers coordinate purely through marker files on the shared filesystem: each
atomically claims an unrun replicate, runs generax in place (GeneRax resumes
from its own checkpoint if the replicate was interrupted), records the resulting
tree, and moves on. There is no shared-memory state, no host discovery, and no
mpirun.

Marker files, per replicate directory:

  running   : the replicate is claimed and being worked on. Its mtime is a
              heartbeat; a `running` marker that has not been touched in
              `_STALE_SECONDS` is assumed orphaned (its crawler died) and may be
              reclaimed by another crawler.
  completed : the replicate finished successfully; `bs-tree.newick` holds its
              gene tree.
  failed    : the replicate failed `_MAX_ATTEMPTS` times and is excluded.
  skipped   : the calculation converged before this replicate ran; it will not
              be run.
  attempts  : integer count of failed attempts so far (drives the retry policy).
"""

import topiary
from topiary._private import Supervisor
from topiary._private.interface import rmtree
from topiary.raxml import check_convergence
from topiary.raxml._raxml import run_raxml

from topiary.generax._generax import GENERAX_BINARY
from topiary.raxml import RAXML_BINARY
from topiary.generax._reconcile_bootstrap import _launch_replicate
from topiary.generax._reconcile_bootstrap import _get_timeout_config
from topiary.generax._reconcile_bootstrap import _compute_replicate_timeout
from topiary.generax._reconcile_bootstrap import _build_replicate_dirs

import os
import glob
import time
import shlex
import shutil
import socket
import tarfile
import threading


# Seconds between heartbeat touches of a claimed replicate's `running` marker.
_HEARTBEAT_INTERVAL = 30.0

# A `running` marker not touched in this many seconds is treated as orphaned
# (its crawler died) and may be reclaimed. Must be comfortably larger than
# _HEARTBEAT_INTERVAL so a slow-but-alive crawler is never declared dead.
_STALE_SECONDS = 300.0

# Total attempts allowed per replicate before it is marked `failed`
# (1 initial attempt + 2 retries).
_MAX_ATTEMPTS = 3

# Seconds a crawler waits before re-scanning when there is nothing claimable but
# work is still in flight on other crawlers.
_POLL_SECONDS = 15.0

# Relative path of the generax result tree within a replicate directory.
_RESULT_TREE = os.path.join("result", "results", "reconcile", "geneTree.newick")

# Coordination marker/lock names.
_SETUP_LOCK = ".topiary-bootstrap-setup.lock"
_CRAWL_READY = ".crawl-ready"
_AGGREGATE_LOCK = ".topiary-bootstrap-aggregate.lock"
_AGGREGATE_DONE = ".aggregate-complete"

# How long a follower waits for setup before giving up.
_SETUP_TIMEOUT = 3600.0


def crawler_id():
    """
    Build an identifier for this crawler process, unique across the allocation.
    """
    return f"{socket.gethostname()}:{os.getpid()}"


def replicate_dirs(replicate_dir):
    """
    Return the sorted list of replicate subdirectories (00001, 00002, ...).

    Parameters
    ----------
    replicate_dir : str
        directory containing the numbered replicate directories

    Returns
    -------
    dirs : list
        sorted list of replicate directory paths
    """

    dirs = [d for d in glob.glob(os.path.join(replicate_dir, "0*"))
            if os.path.isdir(d)]
    dirs.sort()
    return dirs


def _atomic_create(path, content):
    """
    Atomically create `path` (failing if it already exists) and write `content`.

    Parameters
    ----------
    path : str
        file to create
    content : str
        text to write into the file

    Returns
    -------
    created : bool
        True if we created the file; False if it already existed.
    """

    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False

    try:
        os.write(fd, f"{content}\n".encode())
    finally:
        os.close(fd)

    return True


def _is_terminal(d):
    """
    Whether replicate directory `d` is in a terminal state (done, failed, or
    skipped) and therefore not claimable.
    """

    for marker in ("completed", "failed", "skipped"):
        if os.path.exists(os.path.join(d, marker)):
            return True
    return False


def claim_replicate(d, cid, stale_seconds=_STALE_SECONDS):
    """
    Try to claim replicate directory `d` for this crawler.

    A claim is an atomically created `running` marker. If `d` is terminal it is
    not claimable. If it is already claimed but the claim is stale (its heartbeat
    is older than `stale_seconds`), the stale claim is stolen via an atomic
    rename so that at most one crawler wins.

    Parameters
    ----------
    d : str
        replicate directory to claim
    cid : str
        this crawler's id (written into the marker)
    stale_seconds : float
        age past which a claim is considered orphaned

    Returns
    -------
    claimed : bool
        whether this crawler now owns `d`.
    """

    if _is_terminal(d):
        return False

    running = os.path.join(d, "running")

    # Fresh claim.
    if _atomic_create(running, cid):
        # Guard against the tiny race where the replicate was completed between
        # the terminal check and the claim.
        if _is_terminal(d):
            _safe_remove(running)
            return False
        return True

    # Already claimed -- is the claim stale?
    try:
        age = time.time() - os.path.getmtime(running)
    except OSError:
        # Marker vanished (owner finished/released); try once more, fresh.
        return _atomic_create(running, cid) and not _is_terminal(d)

    if age < stale_seconds:
        return False

    # Stale claim: steal it with an atomic rename. Exactly one crawler's rename
    # of this inode can succeed; the rest get FileNotFoundError.
    stolen = f"{running}.steal.{cid.replace(os.sep, '_')}"
    try:
        os.rename(running, stolen)
    except OSError:
        return False

    _safe_remove(stolen)
    if _is_terminal(d):
        return False
    return _atomic_create(running, cid)


def _safe_remove(path):
    """Remove `path`, ignoring the case where it is already gone."""
    try:
        os.remove(path)
    except OSError:
        pass


def _claim_lock(path, cid, stale_seconds=_STALE_SECONDS):
    """
    Atomically claim a coordination lock file, stealing it (via atomic rename) if
    the existing lock's heartbeat is older than `stale_seconds`. Returns True if
    this crawler now holds the lock. Unlike `claim_replicate` there are no
    terminal states -- this is a plain mutual-exclusion lock.
    """

    if _atomic_create(path, cid):
        return True

    try:
        age = time.time() - os.path.getmtime(path)
    except OSError:
        return _atomic_create(path, cid)

    if age < stale_seconds:
        return False

    stolen = f"{path}.steal.{cid.replace(os.sep, '_')}"
    try:
        os.rename(path, stolen)
    except OSError:
        return False

    _safe_remove(stolen)
    return _atomic_create(path, cid)


def find_ready_calc_dir(parent_dir):
    """
    Return the newest ``*reconciled-tree-bootstraps`` calc directory under
    `parent_dir` that has been marked crawl-ready, or None if none exist yet.
    """

    cands = []
    for c in glob.glob(os.path.join(parent_dir, "*reconciled-tree-bootstraps")):
        if not os.path.isdir(c):
            continue
        if not os.path.exists(os.path.join(c, "working", _CRAWL_READY)):
            continue
        try:
            n = int(os.path.basename(c).split("_")[0])
        except ValueError:
            continue
        cands.append((n, c))

    if len(cands) == 0:
        return None

    cands.sort()
    return cands[-1][1]


def elect_setup(parent_dir,
                build_fn,
                cid=None,
                poll=5.0,
                timeout=_SETUP_TIMEOUT,
                stale_seconds=_STALE_SECONDS):
    """
    Ensure bootstrap setup runs exactly once across all crawlers.

    One crawler wins an atomic setup lock in `parent_dir`, runs ``build_fn()``
    (which must build the calc directory and return its path), and drops a
    crawl-ready marker in ``<calc_dir>/working``. The other crawlers wait for a
    crawl-ready calc directory to appear. The setup lock is heartbeated while
    ``build_fn`` runs, so a slow build is never mistaken for a dead leader; if a
    leader really dies mid-build its lock goes stale and another crawler takes
    over.

    Parameters
    ----------
    parent_dir : str
        directory that holds the calc directories (the previous run directory)
    build_fn : callable
        zero-argument function that builds and returns the calc directory path
    cid : str, optional
        crawler id (defaults to hostname:pid)
    poll : float
        seconds between follower re-checks
    timeout : float
        seconds a follower waits before giving up
    stale_seconds : float
        age past which the setup lock is considered orphaned

    Returns
    -------
    calc_dir : str
        the crawl-ready calc directory
    is_leader : bool
        whether this crawler ran setup
    """

    if cid is None:
        cid = crawler_id()

    # Already set up (e.g. a restart)? Just use it.
    existing = find_ready_calc_dir(parent_dir)
    if existing is not None:
        return existing, False

    lock = os.path.join(parent_dir, _SETUP_LOCK)

    def _try_lead():
        if not _claim_lock(lock, cid, stale_seconds):
            return None
        stop = threading.Event()
        hb = threading.Thread(target=_heartbeat_loop,
                              args=(lock, stop, _HEARTBEAT_INTERVAL),
                              daemon=True)
        hb.start()
        try:
            calc_dir = build_fn()
            _atomic_create(os.path.join(calc_dir, "working", _CRAWL_READY), cid)
            return calc_dir
        except Exception:
            _safe_remove(lock)
            raise
        finally:
            stop.set()

    led = _try_lead()
    if led is not None:
        return led, True

    # Follower: wait for a crawl-ready calc dir (re-electing if the leader dies).
    waited = 0.0
    while waited < timeout:
        found = find_ready_calc_dir(parent_dir)
        if found is not None:
            return found, False
        led = _try_lead()
        if led is not None:
            return led, True
        time.sleep(poll)
        waited += poll

    raise RuntimeError("timed out waiting for bootstrap setup to complete")


def elect_aggregate(calc_dir, cid=None, stale_seconds=_STALE_SECONDS):
    """
    Try to become the single aggregator for `calc_dir`. Returns True if this
    crawler should run aggregation (and it has not already been completed).
    """

    if cid is None:
        cid = crawler_id()

    working = os.path.join(calc_dir, "working")
    if os.path.exists(os.path.join(working, _AGGREGATE_DONE)):
        return False

    return _claim_lock(os.path.join(working, _AGGREGATE_LOCK), cid, stale_seconds)


def mark_aggregate_done(calc_dir, cid=None):
    """Record that aggregation for `calc_dir` has completed."""
    if cid is None:
        cid = crawler_id()
    _atomic_create(os.path.join(calc_dir, "working", _AGGREGATE_DONE), cid)


def finalize_bootstrap(calc_dir, converge_cutoff, raxml_binary, report_fn, cid=None):
    """
    Run the one-time finalization (aggregate supports + report) as the elected
    aggregator, holding the aggregate lock the whole time.

    The aggregate lock is heartbeated for the full duration so a long aggregation
    (globbing/reading hundreds of replicate directories on a cluster filesystem
    and compressing them can take many minutes) is never mistaken for a dead
    aggregator and stolen by a second crawler -- which would race on deleting the
    replicates directory.

    The steps are idempotent so a crash is recoverable by simply re-running:
    aggregation is skipped if the supports already exist, the report is
    regenerated every time (it is cheap and overwrites), and the terminal
    ``.aggregate-complete`` marker is written only after the report succeeds.

    Parameters
    ----------
    calc_dir : str
        the reconciliation-bootstrap calc directory
    converge_cutoff : float
        bootstrap convergence criterion
    raxml_binary : str
        raxml binary to use for --support
    report_fn : callable
        zero-argument callback that generates the pipeline report. Kept as a
        callback so this module does not import the (heavy) reporting stack.
    cid : str, optional
        crawler id (defaults to hostname:pid)
    """

    if cid is None:
        cid = crawler_id()

    lock = os.path.join(calc_dir, "working", _AGGREGATE_LOCK)
    stop = threading.Event()
    hb = threading.Thread(target=_heartbeat_loop,
                          args=(lock, stop, _HEARTBEAT_INTERVAL),
                          daemon=True)
    hb.start()
    try:

        # Aggregation is expensive and destructive (it consumes the replicates
        # directory), so skip it if a previous attempt already produced the
        # supports. But we must still make sure the calculation is finalized
        # (calc_status == "complete"): a previous attempt may have produced the
        # supports and then been interrupted before finalizing, and the report
        # ignores any calc directory that is not complete.
        supports = os.path.join(calc_dir, "output",
                                "reconciled-tree_supports.newick")
        if not os.path.isfile(supports):
            aggregate_bootstrap(calc_dir, converge_cutoff,
                                raxml_binary=raxml_binary)
        else:
            sv = Supervisor(calc_dir)
            if sv.status != "complete":
                sv.finalize(successful=True, plot_if_success=True)

        # Report generation is cheap and idempotent; always (re)generate it.
        report_fn()

        # Only now, after the report has succeeded, mark the calculation done so
        # that a crash before this point is recoverable by re-running.
        mark_aggregate_done(calc_dir, cid)

    finally:
        stop.set()


def _heartbeat_loop(running_path, stop_event, interval):
    """
    Touch `running_path` every `interval` seconds until `stop_event` is set or
    the marker disappears. Runs on a daemon thread while generax executes so
    other crawlers can tell this claim is alive.
    """

    while not stop_event.wait(interval):
        try:
            os.utime(running_path, None)
        except OSError:
            return


def _bump_attempts(d):
    """
    Increment (and return) the failed-attempt counter for replicate `d`. Only the
    crawler that currently owns `d` calls this, so no locking is needed.
    """

    attempts_file = os.path.join(d, "attempts")
    try:
        with open(attempts_file) as f:
            n = int(f.read().strip() or "0")
    except (OSError, ValueError):
        n = 0

    n += 1
    with open(attempts_file, "w") as f:
        f.write(f"{n}\n")

    return n


def run_replicate(d, generax_launch, config, durations, sample_threshold):
    """
    Run generax for a single (already-claimed) replicate directory and record
    the outcome. We are launched with the `running` marker in place and a
    heartbeat keeping it fresh.

    generax is run in place: if the replicate was interrupted earlier, GeneRax
    resumes from the checkpoint files already in its `result/` directory.

    Parameters
    ----------
    d : str
        replicate directory (this crawler owns it)
    generax_launch : str
        launcher prefix prepended to the generax command (e.g. "mpirun -np 8").
    config : dict
        timeout configuration (see `_reconcile_bootstrap._DEFAULT_TIMEOUT_CONFIG`)
    durations : list
        this crawler's record of successful replicate runtimes (for adaptive
        timeouts); appended to on success.
    sample_threshold : int
        number of local successes before switching to the adaptive timeout.

    Returns
    -------
    outcome : str
        one of "completed", "retry", or "failed".
    """

    running = os.path.join(d, "running")

    # Build the generax command: the bare `generax ...` line from run_generax.sh
    # with the (optional) user launcher prefix prepended.
    with open(os.path.join(d, "run_generax.sh")) as f:
        bare = f.read().split("\n")[0].strip().split("&>")[0].split()
    cmd = shlex.split(generax_launch) + bare

    timeout = _compute_replicate_timeout(list(durations), sample_threshold, config)

    # Keep the claim alive while generax runs.
    stop = threading.Event()
    hb = threading.Thread(target=_heartbeat_loop,
                          args=(running, stop, _HEARTBEAT_INTERVAL),
                          daemon=True)
    hb.start()

    failure = None
    elapsed = None
    try:
        start = time.time()
        returncode, timed_out = _launch_replicate(cmd,
                                                   timeout=timeout,
                                                   stdout_path=os.path.join(d, "stdout.log"),
                                                   stderr_path=os.path.join(d, "stderr.log"),
                                                   cwd=d)
        elapsed = time.time() - start

        if timed_out:
            failure = f"exceeded timeout of {timeout:.0f} s and was killed"
        elif not os.path.isfile(os.path.join(d, _RESULT_TREE)):
            failure = f"generax exited with code {returncode} and produced no result tree"
    except Exception as e:
        failure = f"unexpected error: {e}"
    finally:
        stop.set()

    if failure is None:

        # Success: publish the tree, then mark completed and drop the claim.
        _publish_tree(d)
        _atomic_create(os.path.join(d, "completed"), crawler_id())
        _safe_remove(running)
        durations.append(elapsed)
        return "completed"

    # Failure: record it, and either release for retry or give up.
    with open(os.path.join(d, "stderr.log"), "a") as f:
        f.write(f"\ntopiary: replicate {failure}\n")

    n = _bump_attempts(d)
    _safe_remove(running)
    if n >= _MAX_ATTEMPTS:
        _atomic_create(os.path.join(d, "failed"), crawler_id())
        return "failed"

    return "retry"


def _publish_tree(d):
    """
    Copy the generax result tree to `bs-tree.newick` in the replicate directory,
    written atomically (temp file + rename) so a reader never sees a partial
    file.
    """

    src = os.path.join(d, _RESULT_TREE)
    dst = os.path.join(d, "bs-tree.newick")
    tmp = f"{dst}.tmp.{os.getpid()}"
    shutil.copy(src, tmp)
    os.replace(tmp, dst)


def assemble_bs_trees(replicate_dir, out_path=None):
    """
    Gather the per-replicate `bs-tree.newick` files from every completed
    replicate into a single multi-tree newick.

    Parameters
    ----------
    replicate_dir : str
        directory containing the replicate directories
    out_path : str, optional
        if given, write the combined trees here.

    Returns
    -------
    trees : list
        list of newick strings, one per completed replicate.
    """

    trees = []
    for d in replicate_dirs(replicate_dir):
        if not os.path.exists(os.path.join(d, "completed")):
            continue
        tree_file = os.path.join(d, "bs-tree.newick")
        if not os.path.isfile(tree_file):
            continue
        with open(tree_file) as f:
            t = f.read().strip()
        if t != "":
            trees.append(t)

    if out_path is not None:
        with open(out_path, "w") as f:
            for t in trees:
                f.write(f"{t}\n")

    return trees


def check_bootstrap_convergence(replicate_dir, converge_cutoff, workdir=None):
    """
    Assemble the completed replicate trees and test bootstrap convergence.

    Parameters
    ----------
    replicate_dir : str
        directory containing the replicate directories
    converge_cutoff : float
        bootstrap convergence criterion (passed to --bs-cutoff)
    workdir : str, optional
        directory in which to write the temporary combined newick. Defaults to
        `replicate_dir`.

    Returns
    -------
    converged : bool
        whether the bootstraps have converged
    df : pandas.DataFrame or None
        convergence report, or None if there were too few trees to test.
    """

    if workdir is None:
        workdir = replicate_dir

    trees = assemble_bs_trees(replicate_dir)
    if len(trees) < 2:
        return False, None

    tmp = os.path.join(workdir, f"convergence-check.{os.getpid()}.newick")
    with open(tmp, "w") as f:
        for t in trees:
            f.write(f"{t}\n")

    try:
        converged, df = check_convergence(tmp, converge_cutoff=converge_cutoff)
    finally:
        _safe_remove(tmp)

    return converged, df


def mark_remaining_skipped(replicate_dir):
    """
    Once converged, drop a `skipped` marker on every replicate that has not been
    claimed or finished so no crawler starts new work.
    """

    for d in replicate_dirs(replicate_dir):
        if _is_terminal(d):
            continue
        if os.path.exists(os.path.join(d, "running")):
            continue
        _atomic_create(os.path.join(d, "skipped"), "converged")


def all_terminal(replicate_dir):
    """
    Whether every replicate is in a terminal state (completed/failed/skipped).
    """

    dirs = replicate_dirs(replicate_dir)
    if len(dirs) == 0:
        return False
    return all(_is_terminal(d) for d in dirs)


def setup_bootstrap(input_bootstrap_dir,
                    calc_dir,
                    converge_cutoff,
                    allow_horizontal_transfer=None,
                    generax_binary=GENERAX_BINARY):
    """
    Build the reconciliation-bootstrap calc directory and its replicate
    directories from a completed ml_bootstrap run. This is the one-time setup
    step; it is run by a single crawler (the setup leader elected by
    `elect_setup`).

    Parameters
    ----------
    input_bootstrap_dir : str
        the previous ml_bootstrap calc directory (holds
        output/bootstrap_replicates)
    calc_dir : str
        name of the reconciliation-bootstrap calc directory to create
    converge_cutoff : float
        bootstrap convergence criterion (stored in run parameters)
    allow_horizontal_transfer : bool, optional
        whether to allow horizontal transfer; if None use whatever the previous
        run recorded, defaulting to True.
    generax_binary : str
        generax binary to use when building the replicate directories

    Returns
    -------
    calc_dir : str
        the calc directory that was built
    """

    supervisor = Supervisor(calc_dir=input_bootstrap_dir)
    supervisor.create_calc_dir(calc_dir=calc_dir, calc_type="reconcile_bootstrap")

    # Must be started from an ml_bootstrap run.
    try:
        prev_calc_type = supervisor.previous_entries[-1]["calc_type"]
    except (IndexError, KeyError, TypeError):
        prev_calc_type = None
    if prev_calc_type != "ml_bootstrap":
        err = "\nbootstrap reconciliation can only be started from a previous\n"
        err += "'ml_bootstrap' run.\n\n"
        raise ValueError(err)

    prev_calc_dir = supervisor.get_previous_calc_dir(-1)
    bs_dir = os.path.join(prev_calc_dir, "output", "bootstrap_replicates")
    if not os.path.isdir(bs_dir):
        err = f"\ninput directory '{prev_calc_dir}' does not have an\n"
        err += "output/bootstrap_replicates directory.\n\n"
        raise FileNotFoundError(err)

    if allow_horizontal_transfer is None:
        if "allow_horizontal_transfer" not in supervisor.run_parameters:
            supervisor.run_parameters["allow_horizontal_transfer"] = True
    else:
        supervisor.update("allow_horizontal_transfer", bool(allow_horizontal_transfer))

    supervisor.check_required(required_values=["model", "allow_horizontal_transfer"],
                              required_files=["alignment.phy", "dataframe.csv",
                                              "gene-tree.newick",
                                              "reconciled-tree.newick"])

    if supervisor.species_tree is None:
        species_tree, _ = topiary.df_to_species_tree(supervisor.df)
        species_tree_out = os.path.join(supervisor.input_dir, "species-tree.newick")
        species_tree.write(outfile=species_tree_out, parser=5)
        supervisor.update("species_tree", species_tree_out)

    supervisor.update("converge_cutoff", converge_cutoff)
    allow_ht = supervisor.run_parameters["allow_horizontal_transfer"]

    starting = os.getcwd()
    os.chdir(supervisor.working_dir)
    try:
        _build_replicate_dirs(df=supervisor.df,
                              model=supervisor.model,
                              gene_tree=supervisor.gene_tree,
                              species_tree=supervisor.species_tree,
                              allow_horizontal_transfer=allow_ht,
                              seed=supervisor.seed,
                              bootstrap_directory=bs_dir,
                              overwrite=False,
                              generax_binary=generax_binary)
    finally:
        os.chdir(starting)

    return calc_dir


def aggregate_bootstrap(calc_dir, converge_cutoff, raxml_binary=RAXML_BINARY):
    """
    Combine the completed bootstrap replicates into branch supports on the
    reconciled tree. This is the one-time aggregation step, run by a single
    crawler (elected by `elect_aggregate`) once every replicate is terminal.

    Parameters
    ----------
    calc_dir : str
        the reconciliation-bootstrap calc directory
    converge_cutoff : float
        bootstrap convergence criterion
    raxml_binary : str
        raxml binary to use for --support

    Returns
    -------
    plot : toyplot.canvas or None
    """

    supervisor = Supervisor(calc_dir)
    starting = os.getcwd()
    os.chdir(supervisor.working_dir)

    replicate_dir = "replicates"

    # Convergence report from the completed replicates.
    converged, df = check_bootstrap_convergence(replicate_dir, converge_cutoff)
    if df is not None:
        df.to_csv("bootstrap-convergence-report.csv")
        supervisor.stash("bootstrap-convergence-report.csv")
    supervisor.update("bootstrap_converged", bool(converged))

    # Assemble the per-replicate trees and map them onto the reconciled tree as
    # branch supports.
    bs_trees = os.path.join(replicate_dir, "bs-trees.newick")
    n_trees = len(assemble_bs_trees(replicate_dir, out_path=bs_trees))
    supervisor.event("Combining bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     num_trees=n_trees)

    run_raxml(run_directory="combine_with_raxml",
              algorithm="--support",
              tree_file=supervisor.reconciled_tree,
              num_threads=1,
              log_to_stdout=False,
              suppress_output=True,
              other_files=[bs_trees],
              other_args=["--bs-trees", "bs-trees.newick", "--redo"],
              raxml_binary=raxml_binary)

    supervisor.stash(os.path.join("combine_with_raxml", "tree.newick.raxml.support"),
                     target_name="reconciled-tree_supports.newick")
    supervisor.stash(os.path.join(replicate_dir, "00001", "species_tree.newick"),
                     "species-tree.newick")

    # Compress the (large) replicates directory and drop it. This is best-effort
    # cleanup: on a big run over a network filesystem it can be slow or get
    # interrupted, and it must never prevent the calculation from being finalized
    # below. Guarded so a re-run after the replicates were already consumed is a
    # no-op.
    if os.path.isdir("replicates"):
        try:
            print("\nCompressing replicates.\n", flush=True)
            f = tarfile.open("replicates.tar.gz", "w:gz")
            f.add("replicates")
            f.close()
            rmtree("replicates")
        except Exception as e:
            print(f"\nWARNING: could not compress the replicates directory "
                  f"({e}). Leaving it in place; the supports have already been "
                  f"computed.\n", flush=True)

    msg = "For more information on the reconciliation events (orthgroups,\n"
    msg += "event counts, full nhx files, etc.) please check the maximum\n"
    msg += "likelihood reconciliation output directory that was used as\n"
    msg += "input for this bootstrap calculation.\n"
    with open(os.path.join(supervisor.output_dir, "reconciliations.txt"), "w") as fh:
        fh.write(msg)

    os.chdir(starting)
    return supervisor.finalize(successful=True, plot_if_success=True)


def crawl(replicate_dir,
          generax_launch="",
          converge_cutoff=0.03,
          timeout_config=None,
          cid=None):
    """
    Run the crawler loop: repeatedly claim and run replicates until none remain
    to be done. Safe to run as many independent copies as desired (across nodes)
    pointing at the same `replicate_dir`.

    Parameters
    ----------
    replicate_dir : str
        directory containing the replicate directories
    generax_launch : str, default=""
        launcher prefix prepended to each generax command.
    converge_cutoff : float, default=0.03
        bootstrap convergence criterion.
    timeout_config : dict, optional
        overrides for the per-replicate timeout configuration.
    cid : str, optional
        crawler id (defaults to hostname:pid).

    Returns
    -------
    n_run : int
        number of replicates this crawler ran to completion.
    """

    if cid is None:
        cid = crawler_id()

    config = _get_timeout_config(timeout_config)
    durations = []
    n_run = 0

    while True:

        # Grab the first claimable replicate.
        claimed = None
        for d in replicate_dirs(replicate_dir):
            if claim_replicate(d, cid):
                claimed = d
                break

        if claimed is None:

            # We are done if every replicate is terminal, or if the replicates
            # directory has already been consumed by the aggregator (it compresses
            # and deletes it). Otherwise other crawlers are still working -- wait
            # and re-scan (a stale claim may free up).
            if not os.path.isdir(replicate_dir) or all_terminal(replicate_dir):
                break
            time.sleep(_POLL_SECONDS)
            continue

        outcome = run_replicate(claimed,
                                generax_launch=generax_launch,
                                config=config,
                                durations=durations,
                                sample_threshold=1)
        if outcome == "completed":
            n_run += 1

            # A completion may have tipped us over the convergence threshold;
            # if so, short-circuit the remaining replicates.
            converged, _ = check_bootstrap_convergence(replicate_dir,
                                                        converge_cutoff)
            if converged:
                mark_remaining_skipped(replicate_dir)

    return n_run
