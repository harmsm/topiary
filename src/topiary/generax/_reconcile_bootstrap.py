"""
Low-level helpers shared by the reconciliation-bootstrap crawler: building the
per-replicate generax directories and running a single generax replicate with a
timeout. The multi-replicate orchestration lives in `topiary.generax._crawl`.
"""

import topiary
from topiary.generax._generax import setup_generax
from topiary.generax._generax import run_generax

from tqdm.auto import tqdm

import os
import glob
import shutil
import subprocess


# Default per-replicate timeout configuration (all times in seconds).
#
# factor  : a replicate is killed if it runs longer than factor * (longest
#           replicate this crawler has seen so far).
# ceiling : timeout to use before we have a runtime estimate (the first
#           replicate a crawler runs).
# floor   : minimum timeout, so fast replicates are not killed by filesystem /
#           scheduler / launcher-startup jitter.
_DEFAULT_TIMEOUT_CONFIG = {"factor": 3.0,
                           "ceiling": 24 * 60 * 60.0,
                           "floor": 300.0}


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


def _launch_replicate(cmd,timeout,stdout_path,stderr_path,cwd=None):
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
    cwd : str, optional
        working directory in which to run the command. The command's relative
        paths (control.txt, --prefix result, ...) resolve here. stdout_path and
        stderr_path are opened relative to the caller's directory, not `cwd`.

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
                                cwd=cwd,
                                env=os.environ.copy())
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
                          log_to_stdout=False,
                          suppress_output=True,
                          write_to_script="run_generax.sh")

    os.chdir(starting_dir)

    return replicate_dir
