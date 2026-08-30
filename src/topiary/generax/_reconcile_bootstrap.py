"""
Reconcile gene and species trees using generax with bootstrap replicates
of the gene tree and alignments.
"""

import topiary

from topiary._private import Supervisor
from topiary._private import threads
from topiary._private import mpi
from topiary._private import run_cleanly
from topiary._private.interface import rmtree
from topiary._private.mpi import get_mpi_env

from topiary.raxml._raxml import run_raxml
from topiary.raxml import check_convergence
from topiary.generax._generax import setup_generax
from topiary.generax._generax import run_generax
from topiary.generax._generax import GENERAX_BINARY


from tqdm.auto import tqdm

import os
import sys
import glob
import shutil
import copy
import tarfile
import random
import string
import subprocess
import time
import pathlib
import multiprocessing as mp


# Default configuration for per-replicate timeouts and the failure circuit
# breaker. All times are in seconds. These can be overridden by the caller via
# the ``timeout_config`` argument threaded down from the pipeline entry point.
#
# factor            : a replicate is killed if it runs longer than
#                     factor * (longest replicate seen so far).
# ceiling           : timeout to use before we have enough completed replicates
#                     to estimate a runtime (i.e. for the very first block of
#                     replicates). Also the maximum a first-block replicate is
#                     allowed to run before we give up on the whole calculation.
# floor             : minimum timeout, so that fast replicates are not killed by
#                     filesystem/scheduler/MPI-startup jitter on a busy cluster.
# max_failed_fraction : abort the whole calculation if more than this fraction
#                     of replicates fail (after `max_failed_floor` failures).
# max_failed_floor  : never trip the circuit breaker until at least this many
#                     replicates have failed (keeps small runs from aborting on
#                     one or two transient failures).
_DEFAULT_TIMEOUT_CONFIG = {"factor":3.0,
                           "ceiling":24*60*60.0,
                           "floor":300.0,
                           "max_failed_fraction":0.1,
                           "max_failed_floor":5}

# How long (seconds) to wait on the shared multiprocessing lock before treating
# it as orphaned (e.g. a sibling worker was hard-killed while holding it). This
# is generously long relative to the sub-second file operations the lock
# protects, so honest contention never trips it.
_LOCK_TIMEOUT = 300.0


class _LocalValue:
    """
    Minimal stand-in for a ``multiprocessing.Manager().Value`` with a ``.value``
    attribute. Used when the calculation runs single-threaded (in-process) and a
    real shared value would be needless overhead.
    """
    def __init__(self,value=0):
        self.value = value


def _get_timeout_config(timeout_config):
    """
    Merge a (possibly partial or None) timeout_config dict on top of the
    defaults.

    Parameters
    ----------
    timeout_config : dict or None
        user-supplied overrides for `_DEFAULT_TIMEOUT_CONFIG`

    Returns
    -------
    config : dict
        complete configuration dictionary
    """

    config = dict(_DEFAULT_TIMEOUT_CONFIG)
    if timeout_config is not None:
        for k in timeout_config:
            if k not in _DEFAULT_TIMEOUT_CONFIG:
                err = f"\nunrecognized timeout_config key '{k}'. allowed keys:\n"
                err += f"{list(_DEFAULT_TIMEOUT_CONFIG.keys())}\n\n"
                raise ValueError(err)
            config[k] = timeout_config[k]

    return config


def _compute_replicate_timeout(durations,
                               sample_threshold,
                               config):
    """
    Decide how long to allow a single generax replicate to run before killing
    it.

    Parameters
    ----------
    durations : list
        wall-clock run times (seconds) of replicates that have completed
        successfully so far.
    sample_threshold : int or None
        number of completed replicates required before we trust the observed
        runtimes and switch from the (long) ceiling to an adaptive timeout. If
        None, always use the ceiling.
    config : dict
        timeout configuration (see `_DEFAULT_TIMEOUT_CONFIG`).

    Returns
    -------
    timeout : float
        number of seconds to allow the replicate to run.
    """

    # Not enough data yet: fall back to the (long) ceiling. This covers the very
    # first block of replicates, where we have no runtime estimate.
    if sample_threshold is None or len(durations) < sample_threshold:
        return config["ceiling"]

    # Adaptive: some multiple of the slowest replicate observed so far, but never
    # less than the floor (protects fast replicates from cluster jitter).
    return max(config["floor"], config["factor"] * max(durations))


def _should_abort(num_failed,
                  num_total,
                  config):
    """
    Decide whether so many replicates have failed that the whole calculation
    should be aborted (rather than silently producing supports from a broken
    run).

    Parameters
    ----------
    num_failed : int
        number of replicates that have failed so far.
    num_total : int
        total number of replicates.
    config : dict
        timeout configuration (see `_DEFAULT_TIMEOUT_CONFIG`).

    Returns
    -------
    abort : bool
        whether or not to abort.
    """

    # Never abort until we have accumulated a meaningful number of failures.
    if num_failed < config["max_failed_floor"]:
        return False

    if num_total is None or num_total <= 0:
        return False

    return num_failed > config["max_failed_fraction"] * num_total


def _terminate_process(proc):
    """
    Stop a launched replicate subprocess (an ``mpirun`` invocation).

    We signal ``mpirun`` itself -- SIGTERM first, so OpenMPI can forward it to
    its ranks and shut them down cleanly, then SIGKILL if it does not exit
    promptly. We deliberately do NOT kill a whole process group: doing so would
    require launching ``mpirun`` in its own session (``start_new_session=True``),
    which interferes with OpenMPI's SLURM launcher (``--mca plm slurm``) and
    serializes what should be parallel replicates. Relying on mpirun's own signal
    forwarding keeps replicates running concurrently at the small cost of a rare
    orphaned rank if mpirun cannot clean up.

    Parameters
    ----------
    proc : subprocess.Popen
        process to terminate.
    """

    # Ask mpirun to shut down (it forwards the signal to its ranks).
    try:
        proc.terminate()
    except (ProcessLookupError,OSError):
        return

    # Give it a chance to exit cleanly.
    try:
        proc.wait(timeout=30)
        return
    except subprocess.TimeoutExpired:
        pass

    # Escalate.
    try:
        proc.kill()
    except (ProcessLookupError,OSError):
        pass


def _launch_replicate(cmd,timeout,stdout_path,stderr_path):
    """
    Run a single replicate command as a subprocess, capturing stdout/stderr to
    files and enforcing a timeout.

    stdout/stderr are redirected to files rather than captured via pipes. This
    avoids the classic MPI deadlock in which mpirun exits but an orphaned rank
    keeps the stdout/stderr pipe open, leaving ``subprocess`` blocked forever
    waiting for EOF. With file redirection the call only waits on the direct
    child (mpirun) exiting.

    The subprocess is launched in the same session/process group as this worker
    (no ``start_new_session``): detaching mpirun into its own session breaks
    OpenMPI's SLURM launcher and serializes otherwise-parallel replicates. On
    timeout we therefore signal mpirun directly rather than killing a group.

    Parameters
    ----------
    cmd : list
        subprocess-style command to run.
    timeout : float
        number of seconds to allow the command to run before killing it.
    stdout_path : str
        file to which stdout is written.
    stderr_path : str
        file to which stderr is written.

    Returns
    -------
    returncode : int or None
        return code of the process (None if it could not be reaped).
    timed_out : bool
        whether or not the process was killed because it hit the timeout.
    """

    stdout_f = open(stdout_path,"w")
    stderr_f = open(stderr_path,"w")
    try:
        proc = subprocess.Popen(cmd,
                                stdout=stdout_f,
                                stderr=stderr_f,
                                env=mpi.get_mpi_env())
    finally:
        # The child has its own duplicated file descriptors; we can close ours.
        stdout_f.close()
        stderr_f.close()

    try:
        proc.wait(timeout=timeout)
        return proc.returncode, False

    except subprocess.TimeoutExpired:

        # Signal mpirun to shut down (it forwards to its ranks), then reap.
        _terminate_process(proc)
        try:
            proc.wait(timeout=30)
        except subprocess.TimeoutExpired:
            pass

        return proc.returncode, True


def _progress_bar(replicate_dir):
    """
    Check to see how far along the calculation is.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    """

    total_calcs = len(glob.glob(os.path.join(replicate_dir,"0*")))

    num_complete = 0
    skipped_added = False
    with tqdm(total=total_calcs) as pbar:
        while num_complete < total_calcs:

            num_complete = len(glob.glob(os.path.join(replicate_dir,"0*","completed")))
            num_running = len(glob.glob(os.path.join(replicate_dir,"0*","running")))
            num_skipped = len(glob.glob(os.path.join(replicate_dir,"0*","skipped")))

            # Num skipped just added in
            if num_skipped > 0 and not skipped_added:

                skipped_added = True
                total_calcs = num_complete + num_running

                pbar.reset(total_calcs)
                pbar.update(num_complete)
                pbar.refresh()

            pbar.n = num_complete
            pbar.refresh()
            time.sleep(0.5)

def _check_convergence(replicate_dir,
                       converge_cutoff,
                       lock=None):
    """
    Check for convergence and indicate to all threads to terminate if true.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    lock : multiprocessing.Manager().Lock(), optional
        lock to control file

    Returns
    -------
    converged : bool
        whether or not the bootstraps are converged
    df : pandas.DataFrame
        dataframe holding convergence results
    """

    if lock is None:
        lock = threads.MockLock()

    # Grab copy of the newick file with lock to avoid collisons with threads
    # that might be writing to it. Use a bounded wait so an orphaned lock (e.g.
    # a sibling worker hard-killed while holding it) cannot wedge the manager
    # thread forever; if we cannot get the lock, simply report not-converged and
    # try again after the next replicate finishes.
    if not lock.acquire(timeout=_LOCK_TIMEOUT):
        return False, None
    try:
        shutil.copy(os.path.join(replicate_dir,"bs-trees.newick"),"tmp.newick")
    finally:
        lock.release()

    # Don't do calculation of the number of trees is only one
    f = open('tmp.newick')
    lines = [line for line in f.readlines() if line.strip() != ""]
    f.close()

    if len(lines) < 2:
        os.remove("tmp.newick")
        return False, None

    # Check for convergence
    converged, df = check_convergence("tmp.newick",
                                      converge_cutoff=converge_cutoff)

    # Delete temporary newick file
    os.remove("tmp.newick")

    # If calculation has converged...
    if converged:

        # Get list of directories
        dirs = [os.path.join(replicate_dir,d) for d in os.listdir(replicate_dir)]
        dirs = [d for d in dirs if os.path.isdir(d)]
        dirs.sort()

        num_completed = 0
        num_running = 0

        # Lock to prevent other threads from starting a calculation or updating
        # the progress bar while this is running. Use a bounded wait so an
        # orphaned lock cannot wedge the manager thread; if we cannot get it,
        # skip marking (the remaining workers will still terminate naturally as
        # they exhaust the directories).
        got_lock = lock.acquire(timeout=_LOCK_TIMEOUT)
        if got_lock:
            try:

                # Go through all directories
                for d in dirs:

                    # Is calculation done?
                    if os.path.isfile(os.path.join(d,"completed")):
                        num_completed += 1

                    # if not...
                    else:

                        # If calculation running?
                        if os.path.isfile(os.path.join(d,"running")):
                            num_running += 1

                        # if not, write file telling other threads to skip
                        else:
                            pathlib.Path(os.path.join(d,"skipped")).touch()

            # Release lock
            finally:
                lock.release()

    return converged, df


def _generax_thread_function(replicate_dir,
                             converge_cutoff,
                             is_manager,
                             hosts,
                             lock=None,
                             durations=None,
                             fail_count=None,
                             total_replicates=None,
                             sample_threshold=None,
                             timeout_config=None):
    """
    Run a generax calculation in parallel, checking for and avoiding collisions
    with other workers.

    Parameters
    ----------

    replicate_dir : str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    is_manager : bool
        whether or not this is the manager thread that will check for
        convergence.
    hosts : list
        list of hosts on which to run the calculation. passed to mpirun via
        --hosts ",".join(hosts)
    lock : multiprocessing.Manager().Lock()
        lock to allow multiple threads to access files
    durations : list or multiprocessing.Manager().list(), optional
        shared list to which the wall-clock runtime of each successful replicate
        is appended. Used to set adaptive timeouts. If None, a private list is
        used (no cross-worker sharing).
    fail_count : multiprocessing.Manager().Value or object with a `.value`, optional
        shared counter of failed replicates, used for the failure circuit
        breaker. If None, a private counter is used.
    total_replicates : int, optional
        total number of replicates in the calculation. Used by the failure
        circuit breaker. If None, the circuit breaker never trips.
    sample_threshold : int, optional
        number of successful replicates required before switching from the
        ceiling timeout to an adaptive timeout. If None, always use the ceiling.
    timeout_config : dict, optional
        overrides for `_DEFAULT_TIMEOUT_CONFIG`.
    """

    if lock is None:
        lock = threads.MockLock()
    if durations is None:
        durations = []
    if fail_count is None:
        fail_count = _LocalValue()

    config = _get_timeout_config(timeout_config)

    # Change into replicate_directory and get sorted list of bootstrap
    # replicates.
    os.chdir(replicate_dir)
    dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
    dirs.sort()

    # Construct a base mpirun command that the generax commands will be
    # appended to. If we are running purely on the local node, omit the --host 
    # flag entirely to avoid OpenMPI attempting to SSH to localhost.
    base_cmd = ["mpirun"]
    base_cmd.extend(mpi.get_mpi_flags())
    if mpi._get_mpi_oversubscribe():
        base_cmd.append("--oversubscribe")
    
    base_cmd.extend(["-np", str(len(hosts))])

    if not all([h == "localhost" for h in hosts]):
        base_cmd.extend(["--host",",".join(hosts)])

    # Path to result tree within each directory
    result_tree = os.path.join("result","results","reconcile","geneTree.newick")

    # Go through each directory
    converged = None
    df = None
    for d in dirs:

        # figure out what directory to go into. Use the lock to make sure only
        # one process is staking a claim at a time.
        lock.acquire()
        try:

            # If this directory has already been run
            if os.path.isfile(os.path.join(d,"completed")):
                continue

            # If this directory is already running
            if os.path.isfile(os.path.join(d,"running")):
                continue

            # If this directory has already been set to skip
            if os.path.isfile(os.path.join(d,"skipped")):
                continue
            
            # If this directory failed previously
            if os.path.isfile(os.path.join(d,"failed")):
                continue

            # Stake a claim
            os.chdir(d)
            pathlib.Path("running").touch()

        finally:
            lock.release()

        # If we got all the way here, this directory is ours to run in.

        # Read run_generax.sh script and construct an mpirun command with the
        # contents of the scripts
        f = open("run_generax.sh")
        contents = f.read()
        f.close()

        # Split on "\n" and strip redirect if presents
        bash_cmd = contents.split("\n")[0].strip()
        bash_cmd = bash_cmd.split("&>")[0]
        bash_cmd = bash_cmd.split()
        bash_cmd = [c.strip() for c in bash_cmd if c.strip() != ""]
        cmd = base_cmd[:]
        cmd.extend(bash_cmd)

        # Run the replicate and decide whether it succeeded. `failure_reason` is
        # None on success and a human-readable string on any failure. Everything
        # here is wrapped so that no error -- expected or not -- can leave the
        # directory flagged `running` and wedge the calculation.
        failure_reason = None
        elapsed = None
        try:

            # Decide how long to let this replicate run. Snapshot the durations
            # seen so far (list() works for a plain list or a manager list proxy).
            timeout = _compute_replicate_timeout(list(durations),
                                                 sample_threshold,
                                                 config)

            # Launch the replicate. Output is redirected to files (rather than
            # captured via pipes) to avoid the MPI deadlock where mpirun exits but
            # an orphaned rank keeps the stdout/stderr pipe open. The timeout
            # guarantees we can never wedge on a single hung replicate.
            start_time = time.time()
            returncode, timed_out = _launch_replicate(cmd,
                                                      timeout=timeout,
                                                      stdout_path="stdout.log",
                                                      stderr_path="stderr.log")
            elapsed = time.time() - start_time

            # MPI can throw a non-zero exit code even when generax finished fine,
            # so the presence of the result tree is the real arbiter of success.
            if timed_out:
                failure_reason = (f"replicate exceeded its timeout of "
                                  f"{timeout:.0f} s and was killed.")
            elif not os.path.isfile(result_tree):
                failure_reason = (f"generax exited with code {returncode} and "
                                  f"did not produce a result tree.")
            else:

                # Success: read the tree and append it to the shared bs-trees
                # file under the lock. Use a bounded wait so an orphaned lock
                # cannot wedge this worker; if we cannot get it, fail this
                # replicate rather than block forever.
                f = open(result_tree,'r')
                tree = f.read().strip()
                f.close()

                if lock.acquire(timeout=_LOCK_TIMEOUT):
                    try:
                        f = open(os.path.join("..","bs-trees.newick"),"a")
                        f.write(f"{tree}\n")
                        f.close()
                    finally:
                        lock.release()
                else:
                    failure_reason = ("could not acquire the shared lock to "
                                      "record the result (a sibling worker may "
                                      "have died holding it).")

        except Exception as e:
            failure_reason = f"unexpected error running replicate: {e}"

        # Success: record the runtime and flag the directory complete.
        if failure_reason is None:

            durations.append(elapsed)
            pathlib.Path("completed").touch()
            if os.path.exists("running"):
                os.remove("running")
            os.chdir("..")

        # Failure (crash, timeout, missing tree, or lost lock): annotate the
        # directory, drop the claim, and move on. We are still inside `d` here.
        else:

            try:
                with open("stderr.log","a") as f:
                    f.write(f"\ntopiary: {failure_reason}\n")
            except OSError:
                pass

            pathlib.Path("failed").touch()
            if os.path.exists("running"):
                os.remove("running")
            os.chdir("..")

            # Track the failure and, if too many replicates have failed, abort
            # the whole calculation rather than produce garbage supports.
            if lock.acquire(timeout=_LOCK_TIMEOUT):
                try:
                    fail_count.value += 1
                    num_failed = fail_count.value
                finally:
                    lock.release()
            else:
                num_failed = fail_count.value

            w = f"\nWARNING: bootstrap replicate {d} failed ({failure_reason})\n"
            w += f"Check {os.path.abspath(os.path.join(d,'stderr.log'))} for details.\n"
            print(w, flush=True)

            if _should_abort(num_failed,total_replicates,config):
                err = f"\n{num_failed} bootstrap replicates have failed, which\n"
                err += "exceeds the allowed failure fraction "
                err += f"({config['max_failed_fraction']}). Aborting. This usually\n"
                err += "indicates a systemic problem (a bad node, an MPI\n"
                err += "misconfiguration, or an input problem) rather than\n"
                err += "isolated replicate failures.\n\n"
                raise RuntimeError(err)

            continue

        # For the manager thread, check for convergence
        if is_manager:
            converged, df = _check_convergence(replicate_dir=".",
                                               converge_cutoff=converge_cutoff,
                                               lock=lock)
            if converged:
                break

    os.chdir("..")

    if is_manager:
        return (converged, df)

    return None



def _build_replicate_dirs(df,
                          model,
                          gene_tree,
                          species_tree,
                          allow_horizontal_transfer,
                          seed,
                          bootstrap_directory,
                          overwrite,
                          generax_binary):
    """
    Prepare for generax bootstrap by constructing individual directories that
    hold each bootstrap replicate. Do a dummy run of generax in each one,
    creating a file called run_generax.sh that can be invoked later.

    Parameters
    ----------
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str, optional
        gene tree in newick format.
    species_tree : str, optional
        species tree in newick format.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconciliation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap_directory : str
        directory with gene tree bootstrap replicates
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    replicate_dir : str
        directory containing replicate directories
    """

    starting_dir = os.getcwd()

    # dataframe for boostrap replicates --> only need keep == True
    bs_df = df.loc[df.keep,:]

    # Read species tree
    if species_tree is not None:
        species_tree = species_tree
    else:
        species_tree, dropped = topiary.df_to_species_tree(bs_df,strict=True)
        species_tree.resolve_polytomy()

    # Create replicate directory
    replicate_dir = os.path.abspath(os.path.join("replicates"))
    os.mkdir(replicate_dir)

    # Imported here rather than at the top of the module because importing ete4
    # pulls in scipy, which is slow (see also topiary/__init__.py).
    import ete4 as ete

    # Find and sort bootstrap alignment files
    alignment_files = glob.glob(os.path.join(bootstrap_directory,"*.phy"))
    alignment_files = [(int(os.path.split(a)[-1].split("_")[1].split(".")[0]),a)
                       for a in alignment_files]
    alignment_files.sort()
    alignment_files = [a[1] for a in alignment_files]

    # Load bootstrap trees
    tree_files = glob.glob(os.path.join(bootstrap_directory,"*.newick"))[0]
    trees = []
    with open(tree_files) as f:
        for line in f:
            trees.append(ete.Tree(line.strip()))

    # Sanity check on number of trees vs number of alignments
    if len(trees) != len(alignment_files):
        err = "\nNumber of bootstrap trees does not match the number of bootstrap\n"
        err += "alignments.\n\n"
        raise ValueError(err)

    template_species_tree = None
    template_control_file = None
    template_link_file = None
    template_keep_mask = None

    # For every alignment
    for i in tqdm(range(len(alignment_files))):

        rep_number = f"{(i+1):05d}"

        this_df = bs_df.copy()

        out_dir = os.path.abspath(os.path.join(replicate_dir,rep_number))

        # Write to replicate working directory
        keep_mask = setup_generax(this_df,
                                  trees[i],
                                  model,
                                  out_dir,
                                  keep_mask=template_keep_mask,
                                  species_tree=template_species_tree,
                                  mapping_link_file=template_link_file,
                                  control_file=template_control_file)

        # If this is None, we have only done first iteration. Capture those
        # files to pass into next run. (Saves a ton of time to not have to
        # recalculate species tree, traverse trees for mapping, etc.)
        if template_species_tree is None:
            template_species_tree = os.path.abspath(os.path.join(out_dir,"species_tree.newick"))
            template_control_file = os.path.abspath(os.path.join(out_dir,"control.txt"))
            template_link_file = os.path.abspath(os.path.join(out_dir,"mapping.link"))
            template_keep_mask = keep_mask.copy()

        # Copy the alignment file in, overwriting what was generated by
        # setup_generax
        align_file = os.path.join(out_dir,"alignment.phy")
        shutil.copy(alignment_files[i],align_file)

        os.chdir(replicate_dir)

        # This is a dummy run -- exactly as it would be called in a standard
        # reconciliation -- but now we're writing resulting command to
        # run_generax.sh in the run directory.
        cmd = run_generax(run_directory=out_dir,
                          allow_horizontal_transfer=allow_horizontal_transfer,
                          seed=seed,
                          generax_binary=generax_binary,
                          num_threads=1,
                          log_to_stdout=False,
                          suppress_output=True,
                          write_to_script="run_generax.sh")

    os.chdir(starting_dir)

    return replicate_dir

def _clean_replicate_dir(replicate_dir):
    """
    Remove previous run information.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    """

    running = glob.glob(os.path.join(replicate_dir,"*","running"),recursive=True)
    for f in running:
        os.remove(f)

    skipped = glob.glob(os.path.join(replicate_dir,"*","skipped"),recursive=True)
    for f in skipped:
        os.remove(f)


def _construct_args(replicate_dir,
                    converge_cutoff,
                    num_threads,
                    threads_per_rep,
                    durations=None,
                    fail_count=None,
                    total_replicates=None,
                    timeout_config=None):
    """
    Construct a list of arguments to pass to each thread in the pool.

    Parameters
    ----------
    replicate_dir: str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    num_threads : int
        total number of mpi slots to use for the calculation
    threads_per_rep : int
        number of slots to use per replicate.
    durations : list or multiprocessing.Manager().list(), optional
        shared list of successful-replicate runtimes (for adaptive timeouts).
    fail_count : object with a `.value` attribute, optional
        shared counter of failed replicates (for the circuit breaker).
    total_replicates : int, optional
        total number of replicate directories (for the circuit breaker).
    timeout_config : dict, optional
        overrides for `_DEFAULT_TIMEOUT_CONFIG`.

    Returns
    -------
    kwargs_list : list
        list of dictionaries with kwargs to pass for each calculation
    num_threads : int
        number of total calculations to start via thread_manager
    """

    hosts = mpi.get_hosts(num_threads)

    kwargs_list = []
    for i in range(0,len(hosts),threads_per_rep):

        this_hosts = hosts[i:i+threads_per_rep]

        # Check for node locality
        unique_hosts = set(this_hosts)
        if len(unique_hosts) > 1:
            w = f"\nWARNING: A bootstrap replicate calculation spans multiple compute\n"
            w += f"nodes ({', '.join(list(unique_hosts))}). This can lead to very slow\n"
            w += "performance. Consider setting threads_per_replicate such that\n"
            w += "it is a factor of the number of slots per node.\n"
            print(w, flush=True)

        if i == 0:
            is_manager = True
        else:
            is_manager = False

        kwargs_list.append({"replicate_dir":replicate_dir,
                            "is_manager":is_manager,
                            "hosts":this_hosts,
                            "converge_cutoff":converge_cutoff})

    num_workers = len(kwargs_list)

    # Only start trusting observed runtimes (and switch from the ceiling to an
    # adaptive timeout) once a full first block of workers has each completed at
    # least one replicate. Cap by the total number of replicates so tiny runs
    # still eventually adapt.
    if total_replicates is not None:
        sample_threshold = min(num_workers,total_replicates)
    else:
        sample_threshold = num_workers

    # Inject shared state / timeout configuration into every worker's kwargs.
    for kwargs in kwargs_list:
        kwargs["durations"] = durations
        kwargs["fail_count"] = fail_count
        kwargs["total_replicates"] = total_replicates
        kwargs["sample_threshold"] = sample_threshold
        kwargs["timeout_config"] = timeout_config

    return kwargs_list, num_workers


def _run_bootstrap_calculations(replicate_dir,
                                converge_cutoff,
                                num_threads,
                                threads_per_rep,
                                timeout_config=None):
    """
    Run generax in parallel using mpirun for all directories.

    Parameters
    ----------
    replicate_dir : str
        directory with bootstrap replicates for reconciliations
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    num_threads : int
        number of parallel jobs to start
    threads_per_rep : int
        number of threads to use per replicate. only used if bootstrap = True
    timeout_config : dict, optional
        overrides for `_DEFAULT_TIMEOUT_CONFIG` (per-replicate timeout and
        failure circuit breaker).
    """

    print("\nGenerating reconciliation bootstraps.\n",flush=True)

    # Total number of replicate directories (used by the failure circuit breaker)
    total_replicates = len(glob.glob(os.path.join(replicate_dir,"0*")))

    # Shared state for adaptive timeouts and the failure circuit breaker. When
    # single-threaded, the workers run in-process, so plain objects suffice and
    # we avoid spinning up a manager. When multi-threaded, use manager proxies
    # so the separate worker processes share the same state.
    manager = None
    if num_threads == 1:
        durations = []
        fail_count = _LocalValue(0)
    else:
        manager = mp.Manager()
        durations = manager.list()
        fail_count = manager.Value("i",0)

    # This is a status bar that we spawn on its own thread that will spew onto
    # stderr as the calculations are completed in each directory.
    status_bar = mp.Process(target=_progress_bar,args=(replicate_dir,))
    status_bar.start()

    # This TMPDIR insanity is to prevent MPI from choking because temporary
    # directory names can get too long for it's buffer to handle. (Really.)
    try:
        old_tmpdir = os.environ["TMPDIR"]
    except KeyError:
        old_tmpdir = None
    os.environ["TMPDIR"] = "/tmp/"


    kwargs_list, num_threads = _construct_args(replicate_dir,
                                               converge_cutoff,
                                               num_threads,
                                               threads_per_rep,
                                               durations=durations,
                                               fail_count=fail_count,
                                               total_replicates=total_replicates,
                                               timeout_config=timeout_config)

    # Launch calculation.
    try:
        results = threads.thread_manager(kwargs_list,
                                         _generax_thread_function,
                                         num_threads=num_threads,
                                         progress_bar=False,
                                         pass_lock=True)

    # This clean up step makes sure we kill the status bar thread if the
    # calculation crashes.
    except Exception as e:
        status_bar.kill()
        if manager is not None:
            manager.shutdown()
        raise e

    # If we get here, the job is done whether the status bar is or not. Wait
    # two seconds to let the status bar finish if it's going to so the log
    # reflects the current status of the run.
    time.sleep(2)
    status_bar.kill()

    # Revert TMPDIR if necessary
    if old_tmpdir is not None:
        os.environ["TMPDIR"] = old_tmpdir

    # Only the first thread will spit out results: a tuple with convergence
    # status (True or False) and the dataframe corresponding to last convergence
    # test. Get those results. If df is None, convergence test not run yet;
    # run it.
    converged, df = results[0]
    if df is None:
        converged, df = _check_convergence(replicate_dir=replicate_dir,
                                           converge_cutoff=converge_cutoff)

    if manager is not None:
        manager.shutdown()

    return converged, df

@run_cleanly
def reconcile_bootstrap(df,
                        model,
                        gene_tree,
                        species_tree,
                        reconciled_tree,
                        allow_horizontal_transfer,
                        seed,
                        bootstrap_directory,
                        converge_cutoff,
                        restart,
                        overwrite,
                        supervisor,
                        num_threads,
                        threads_per_rep,
                        generax_binary,
                        raxml_binary,
                        timeout_config=None):
    """
    Reconcile gene and species trees using generax with bootstrap replicates
    of the gene tree and alignments.

    Parameters
    ----------
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str, ete4.Tree, dendropy.tree, optional
        gene tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete4 or dendropy tree, it will be written out with leaf
        names and branch lengths; all other data will be dropped.
    species_tree : str, ete4.Tree, dendropy.tree, optional
        species tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete4 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped.
    reconciled_tree : str, ete4.Tree, dendropy.tree, optional
        reconciled tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete4 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconciliation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap: bool, default=False
        whether or not to do bootstrap replicates. if True, prev_calculation must
        point to a raxml ml_bootstrap run
    converge_cutoff : float, default=0.03
        bootstrap convergence criterion. This is RAxML-NG default, passed to
        --bs-cutoff.
    restart : str, optional
        if specified, should point to replicate_dir for restart. restart the job
        from where it stopped. incompatible with overwrite
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    supervisor : Supervisor, optional
        supervisor instance to keep track of inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1 use all available.
    threads_per_rep : int, default=1
        number of threads to use per replicate.
    generax_binary : str, optional
        what generax binary to use
    timeout_config : dict, optional
        overrides for `_DEFAULT_TIMEOUT_CONFIG`, controlling the per-replicate
        timeout (keys "factor", "ceiling", "floor") and the failure circuit
        breaker (keys "max_failed_fraction", "max_failed_floor").

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    os.chdir(supervisor.working_dir)

    if restart is None:

        # Create stack of directories
        supervisor.event("Creating bootstrap directories.",
                         model=model,
                         gene_tree=gene_tree,
                         species_tree=species_tree,
                         allow_horizontal_transfer=allow_horizontal_transfer,
                         seed=seed,
                         bootstrap_directory=bootstrap_directory,
                         overwrite=overwrite,
                         generax_binary=generax_binary)

        replicate_dir = _build_replicate_dirs(df=df,
                                              model=model,
                                              gene_tree=gene_tree,
                                              species_tree=species_tree,
                                              allow_horizontal_transfer=allow_horizontal_transfer,
                                              seed=seed,
                                              bootstrap_directory=bootstrap_directory,
                                              overwrite=overwrite,
                                              generax_binary=generax_binary)

    else:

        # Restart, making sure any leftover claim/skipped files are removed
        replicate_dir = restart
        _clean_replicate_dir(replicate_dir)



    supervisor.event("Running bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     converge_cutoff=converge_cutoff,
                     num_threads=num_threads,
                     threads_per_rep=threads_per_rep)
    converged, df = _run_bootstrap_calculations(replicate_dir,
                                                converge_cutoff,
                                                num_threads,
                                                threads_per_rep,
                                                timeout_config=timeout_config)

    # Write convergence report and whether this converged or not
    df.to_csv("bootstrap-convergence-report.csv")
    supervisor.stash("bootstrap-convergence-report.csv")
    supervisor.update("bootstrap_converged",bool(converged))

    # Combine bootstrap replicates into a set of supports
    supervisor.event("Combining bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     reconciled_tree=reconciled_tree)

    bs_trees = os.path.join(replicate_dir,"bs-trees.newick")
    cmd = run_raxml(run_directory="combine_with_raxml",
                    algorithm="--support",
                    tree_file=reconciled_tree,
                    num_threads=1,
                    log_to_stdout=False,
                    suppress_output=True,
                    other_files=[bs_trees],
                    other_args=["--bs-trees","bs-trees.newick","--redo"])

    # Copy reconciled tree with suppports into output
    supervisor.stash(os.path.join("combine_with_raxml","tree.newick.raxml.support"),
                     target_name="reconciled-tree_supports.newick")

    # Grab species tree before compressing replicates
    supervisor.stash(os.path.join(replicate_dir,"00001","species_tree.newick"),
                     "species-tree.newick")

    # Compress big, complicated replicates directory and delete
    print("\nCompressing replicates.\n",flush=True)
    f = tarfile.open("replicates.tar.gz","w:gz")
    f.add("replicates")
    f.close()
    rmtree("replicates")

    # Write message indicating where to look for further output
    msg = "For more information on the reconciliation events (orthgroups,\n"
    msg += "event counts, full nhx files, etc.) please check the maximum\n"
    msg += "likelihood reconciliation output directory that was used as\n"
    msg += "input for this bootstrap calculation.\n"

    f = open(os.path.join(supervisor.output_dir,"reconciliations.txt"),"w")
    f.write(msg)
    f.close()

    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
