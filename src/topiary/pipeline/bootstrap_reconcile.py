"""
Final step of the pipeline: compute gene/species-tree reconciliation bootstrap
supports.

Each bootstrap replicate is an independent generax reconciliation. This command
is a re-entrant *crawler*: launch as many copies as you like (e.g. one per
compute node via a SLURM job array), all pointing at the same previous run
directory. They coordinate purely through the filesystem -- one builds the
replicate directories, all of them run replicates, and one assembles the final
supports. There is no MPI orchestration; parallelism across generax processes,
if wanted, is left to the user via ``--generax-launch``.
"""

from topiary.raxml import RAXML_BINARY
from topiary.generax import GENERAX_BINARY
from topiary.generax import _crawl
from topiary._private import installed
from topiary._private import software_requirements
from topiary._private import check
from topiary._private import Supervisor
from topiary._private import run_cleanly

import os
import glob


@run_cleanly
def bootstrap_reconcile(previous_run_dir,
                        generax_launch="",
                        converge_cutoff=0.03,
                        replicate_timeout_factor=3.0,
                        replicate_max_hours=24.0,
                        replicate_min_seconds=300.0,
                        raxml_binary=RAXML_BINARY,
                        generax_binary=GENERAX_BINARY):
    """
    Compute bootstrap branch supports for the gene/species tree reconciliation.

    This command is re-entrant and coordinates through the filesystem, so it is
    safe (and intended) to run many copies at once against the same
    ``previous_run_dir`` -- for example a SLURM job array with one task per node.
    The first crawler builds the replicate directories; every crawler runs
    replicates; the last one assembles the supports and writes the report. Re-run
    it to resume an interrupted calculation -- completed replicates are skipped
    and interrupted ones resume from GeneRax's own checkpoint. To start over,
    delete the ``*_reconciled-tree-bootstraps`` directory first.

    Parameters
    ----------
    previous_run_dir : str
        previous pipeline run directory. Must contain a completed
        ``xx_*bootstraps`` (ml_bootstrap) directory to use as input.
    generax_launch : str, default=""
        launcher prefix prepended to every generax command (e.g. "mpirun -np 8").
        Empty runs each replicate as a single process. GeneRax parallelism is
        MPI-only; if you set this, YOU are responsible for the launcher being
        available and for allocating the resources it needs on each node.
    converge_cutoff : float, default=0.03
        bootstrap convergence criterion (RAxML-NG default, passed to --bs-cutoff).
    replicate_timeout_factor : float, default=3.0
        a replicate is killed if it runs longer than this factor times the
        longest replicate this crawler has seen. Guards against a hung generax.
    replicate_max_hours : float, default=24.0
        per-replicate timeout used before a crawler has a runtime estimate.
    replicate_min_seconds : float, default=300.0
        minimum per-replicate timeout, so fast replicates are not killed by
        filesystem/scheduler jitter.
    raxml_binary : str, optional
        raxml binary to use (for the final --support step)
    generax_binary : str, optional
        generax binary to use
    """

    # --------------------------------------------------------------------------
    # Argument checks

    if not os.path.isdir(previous_run_dir):
        err = f"\nprevious_run_dir '{previous_run_dir}' does not exist\n\n"
        raise FileNotFoundError(err)

    converge_cutoff = check.check_float(converge_cutoff, "converge_cutoff",
                                        minimum_allowed=0, maximum_allowed=1)
    replicate_timeout_factor = check.check_float(replicate_timeout_factor,
                                                 "replicate_timeout_factor",
                                                 minimum_allowed=1.0)
    replicate_max_hours = check.check_float(replicate_max_hours,
                                            "replicate_max_hours",
                                            minimum_allowed=0,
                                            minimum_inclusive=False)
    replicate_min_seconds = check.check_float(replicate_min_seconds,
                                              "replicate_min_seconds",
                                              minimum_allowed=0,
                                              minimum_inclusive=False)

    timeout_config = {"factor": replicate_timeout_factor,
                      "ceiling": replicate_max_hours * 60 * 60,
                      "floor": replicate_min_seconds}

    # --------------------------------------------------------------------------
    # Validate software stack (no mpirun -- topiary no longer orchestrates MPI)

    to_validate = [{"program": "raxml-ng",
                    "binary": raxml_binary,
                    "min_version": software_requirements["raxml-ng"],
                    "must_pass": True},
                   {"program": "generax",
                    "binary": generax_binary,
                    "min_version": software_requirements["generax"],
                    "must_pass": True}]
    installed.validate_stack(to_validate)

    # --------------------------------------------------------------------------
    # Locate the completed ml_bootstrap directory to use as input

    starting_dir = os.getcwd()
    os.chdir(previous_run_dir)
    try:
        _run(previous_run_dir=previous_run_dir,
             generax_launch=generax_launch,
             converge_cutoff=converge_cutoff,
             timeout_config=timeout_config,
             generax_binary=generax_binary,
             raxml_binary=raxml_binary)
    finally:
        os.chdir(starting_dir)


def _input_bootstrap_dir(previous_run_dir):
    """
    Find the highest-numbered ``*bootstraps`` directory that is NOT a
    reconciliation-bootstrap directory (i.e. the ml_bootstrap input). Assumes cwd
    is `previous_run_dir`.
    """

    cands = []
    for b in glob.glob("*bootstraps*"):
        if not os.path.isdir(b):
            continue
        if b.endswith("reconciled-tree-bootstraps"):
            continue
        try:
            n = int(b.split("_")[0])
        except ValueError:
            continue
        cands.append((n, b))

    if len(cands) == 0:
        err = f"\nprevious_run_dir '{previous_run_dir}' does not contain an input\n"
        err += "bootstraps directory (a completed ml_bootstrap run). This is\n"
        err += "required as input to a reconciliation bootstrap calculation.\n\n"
        raise FileNotFoundError(err)

    cands.sort()
    return cands[-1]


def _run(previous_run_dir, generax_launch, converge_cutoff,
         timeout_config, generax_binary, raxml_binary):
    """
    Core re-entrant crawler flow (cwd is `previous_run_dir`).
    """

    cid = _crawl.crawler_id()

    input_num, input_dir = _input_bootstrap_dir(previous_run_dir)

    # Make sure the input calculation completed.
    input_supervisor = Supervisor(input_dir)
    if input_supervisor.status != "complete":
        err = f"\ninput '{previous_run_dir}/{input_dir}' has status\n"
        err += f"'{input_supervisor.status}', not 'complete'. The ml_bootstrap\n"
        err += "calculation must finish before computing reconciliation supports.\n\n"
        raise RuntimeError(err)

    # Setup (runs exactly once across all crawlers): build the replicate dirs.
    def _build():
        calc_dir = f"{input_num + 1:02d}_reconciled-tree-bootstraps"
        return _crawl.setup_bootstrap(input_bootstrap_dir=input_dir,
                                      calc_dir=calc_dir,
                                      converge_cutoff=converge_cutoff,
                                      generax_binary=generax_binary)

    # cwd is already previous_run_dir, so coordinate in "." (the calc directories
    # and the setup lock all live directly under it).
    calc_dir, _is_leader = _crawl.elect_setup(".", _build, cid=cid)

    replicate_dir = os.path.join(calc_dir, "working", "replicates")

    # Crawl: claim and run replicates until every replicate is terminal.
    _crawl.crawl(replicate_dir,
                 generax_launch=generax_launch,
                 converge_cutoff=converge_cutoff,
                 timeout_config=timeout_config,
                 cid=cid)

    # Aggregate: exactly one crawler assembles the supports and writes the
    # report. finalize_bootstrap heartbeats the aggregate lock for the whole
    # (potentially long) operation and only marks the calculation done after the
    # report succeeds, so it is safe against a slow filesystem and recoverable by
    # re-running.
    run_dir = os.getcwd()

    def _report():
        # Imported here rather than at module top because report generation pulls
        # in matplotlib and toytree, which are slow to import.
        from topiary.reports import pipeline_report
        pipeline_report(pipeline_directory=run_dir,
                        output_directory=os.path.join(run_dir, "results"),
                        overwrite=True)

    if _crawl.all_terminal(replicate_dir) and _crawl.elect_aggregate(calc_dir, cid):
        _crawl.finalize_bootstrap(calc_dir,
                                  converge_cutoff=converge_cutoff,
                                  raxml_binary=raxml_binary,
                                  report_fn=_report,
                                  cid=cid)
