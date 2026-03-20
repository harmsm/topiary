"""
Final step on the pipeline. Run replicates in embarrassingly parallel fashion
across compute nodes using MPI.
"""

import topiary
from topiary.raxml import RAXML_BINARY
from topiary.generax import GENERAX_BINARY
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary._private import installed
from topiary._private import software_requirements
from topiary._private.mpi import check_mpi_configuration
from topiary._private import check
from topiary._private import Supervisor
from topiary._private import run_cleanly
from topiary._private.interface import rmtree
from topiary.reports import pipeline_report

import os
import datetime
import glob
import shutil

@run_cleanly
def bootstrap_reconcile(previous_run_dir,
                        num_threads=None,
                        threads_per_replicate=None,
                        converge_cutoff=0.03,
                        restart=False,
                        overwrite=False,
                        raxml_binary=RAXML_BINARY,
                        generax_binary=GENERAX_BINARY):
    """
    Perform a bootstrap branch support calculation for the gene/species tree
    reconciliation portion of the analysis.

    previous_run_dir : str
        previous pipeline run directory. Should have a directory named 
        xx_*bootstraps*, where xx is an integer and * are any value. 
    num_threads : int, optional
        total number of threads (slots, in MPI lingo) to use. If None, topiary
        will infer the number of slots from the environment.
    threads_per_replicate : int, optional
        number of threads to use for each bootstrap replicate. To minimize
        the impact of slow cross-node communication, topiary attempts to 
        run each replicate on a single physical node or processor. 
        threads_per_replicate sets the number of slots to use for each 
        replicate. If this is not specified, topiary will choose a number of
        threads based on the number of slots on each compute node. NOTE: if you 
        manually set threads_per_replicate, choose a number that is a factor of
        the number of slots on each compute node to avoid wasting slots. If you
        have 24 slots per node, you could choose 2, 3, 4, 6, 8, 12, or 24.
    converge_cutoff : float, default=0.03
        bootstrap convergence criterion. This is RAxML-NG default, passed 
        to --bs-cutoff.
    restart : bool, default=False
        restart job from where it stopped in output directory. incompatible with
        overwrite
    overwrite : bool, default=False
        whether or not to overwrite existing output. incompatible with restart.
        This will overwrite an existing 05_reconciliation-bootstraps directory,
        not the rest of the pipeline directory.
    raxml_binary : str, optional
        raxml binary to use
    generax_binary : str, optional
        what generax binary to use
    """

    # Make sure pipeline directory is present
    if not os.path.isdir(previous_run_dir):
        err = f"\nprevious_run_dir '{previous_run_dir}' does not exist\n\n"
        raise FileNotFoundError(err)

    # --------------------------------------------------------------------------
    # Check sanity of num_threads and threads_per_replicate
    
    if num_threads is None:
        num_threads = topiary._private.mpi.get_num_slots()

    num_threads = check.check_int(num_threads,
                                  "num_threads",
                                  minimum_allowed=1)

    if threads_per_replicate is None:
        
        # Get hosts to see how many slots per node
        hosts = topiary._private.mpi.get_hosts(num_threads)
        
        # Get number of slots on the first node (assume homogeneous nodes)
        # This counts how many times the first host appears in the list
        slots_per_node = hosts.count(hosts[0])

        # Find factor of slots_per_node closest to 10
        best_factor = 1
        best_diff = 10
        for i in range(1, slots_per_node + 1):
            if slots_per_node % i == 0:
                diff = abs(i - 10)
                if diff < best_diff:
                    best_diff = diff
                    best_factor = i
                elif diff == best_diff:
                    # If tie (e.g. 7 vs 13), pick smaller
                    if i < best_factor:
                        best_factor = i
        
        threads_per_replicate = best_factor

    threads_per_replicate = check.check_int(threads_per_replicate,
                                            "threads_per_replicate",
                                            minimum_allowed=1)

    # --------------------------------------------------------------------------
    # Check sanity of overwrite, restart, and combination

    overwrite = check.check_bool(overwrite,"overwrite")
    restart = check.check_bool(restart,"restart")

    if overwrite and restart:
        err = "overwrite and restart flags are incompatible.\n"
        raise ValueError(err)

    # --------------------------------------------------------------------------
    # Validate software stack required for this pipeline

    to_validate = [{"program":"raxml-ng",
                    "binary":raxml_binary,
                    "min_version":software_requirements["raxml-ng"],
                    "must_pass":True}]

    to_validate.append({"program":"generax",
                        "binary":generax_binary,
                        "min_version":software_requirements["generax"],
                        "must_pass":True})

    to_validate.append({"program":"mpirun",
                        "min_version":software_requirements["mpirun"],
                        "must_pass":True})

    installed.validate_stack(to_validate)

    # --------------------------------------------------------------------------
    # Validate the previous calculation
    
    os.chdir(previous_run_dir)

    try:

        # All bootstrap directories
        bootstrap_dirs = glob.glob("*bootstraps*")
        if len(bootstrap_dirs) == 0:
            raise ValueError
        
        bootstrap_dirs = [(int(b.split("_")[0]),b) for b in bootstrap_dirs]
        bootstrap_dirs.sort()
        
        # Default behavior: find the highest-numbered bootstrap directory and 
        # assume it's the input. 
        bootstrap_directory = bootstrap_dirs[-1][1]
        dir_counter = bootstrap_dirs[-1][0]
        
        # If we are restarting, we need to be more clever. 
        if restart:

            # If the highest-numbered directory is a reconciled-tree-bootstraps 
            # directory, then WE are the restart. 
            if bootstrap_directory.endswith("reconciled-tree-bootstraps"):
                
                # Check for input *before* this directory. 
                if len(bootstrap_dirs) < 2:
                    err = f"previous_run_dir '{previous_run_dir}' only has a\n"
                    err += "reconciliation bootstrap directory. To restart, there\n"
                    err += "must be an input bootstrap directory as well.\n\n"
                    os.chdir("..")
                    raise RuntimeError(err)
                
                # Find the highest-numbered bootstrap directory that is NOT 
                # a reconciled-tree-bootstraps directory. 
                found_input = False
                for i in range(len(bootstrap_dirs)-2,-1,-1):
                    if not bootstrap_dirs[i][1].endswith("reconciled-tree-bootstraps"):
                        bootstrap_directory = bootstrap_dirs[i][1]
                        found_input = True
                        break
                
                if not found_input:
                    err = f"previous_run_dir '{previous_run_dir}' does not have an\n"
                    err += "input bootstrap directory (other than the reconciliation\n"
                    err += "bootstrap directory itself).\n\n"
                    os.chdir("..")
                    raise RuntimeError(err)

                # The directory we are restarting is the original highest-numbered
                # directory.
                calc_dir = bootstrap_dirs[-1][1]

            # If the highest-numbered directory is NOT a reconciled-tree-bootstraps
            # directory, then we can't restart.
            else:
                err = f"previous_run_dir '{previous_run_dir}' does not have a\n"
                err += "reconciliation bootstrap directory to restart. To start a\n"
                err += "new calculation, do not specify --restart.\n\n"
                os.chdir("..")
                raise RuntimeError(err)

        else:
            calc_dir = f"{dir_counter+1:02d}_reconciled-tree-bootstraps"

    except (ValueError,IndexError):
        err = f"previous_run_dir '{previous_run_dir}' does not have any bootstraps\n"
        err += "directory. This directory is necessary as the input to a\n"
        err += "reconciliation bootstrap calculation.\n\n"
        os.chdir("..")
        raise FileNotFoundError(err)

    # Load calculation and make sure it completed
    supervisor = Supervisor(bootstrap_directory)
    if supervisor.status != "complete":
        err = f"{previous_run_dir}/{bootstrap_directory} exists but has status '{supervisor.status}'\n"
        if supervisor.status == "empty":
            err += "It does not appear this calculation has been run.\n\n"
        elif supervisor.status == "running":
            err += "This job is either still running or crashed.\n\n"
        else:
            err += "This job crashed before completing.\n\n"
        os.chdir("..")
        raise RuntimeError(err)

    # Get number of replicates. Make sure user did not request more slots than
    # replicates.
    num_replicates = len(glob.glob(os.path.join(bootstrap_directory,
                                                "output",
                                                "bootstrap_replicates",
                                                "*.phy")))
    if num_threads > num_replicates:
        print(f"\nWARNING: The number of requested threads (slots: {num_threads}) is\n"
              f"greater than the number of bootstrap replicates ({num_replicates}).\n"
              f"Dropping the number of threads to {num_replicates} to match.\n", flush=True)
        num_threads = num_replicates

    # Now that we've potentially dropped the number of threads to match the
    # number of replicates, check to see if mpirun can actually grab them.
    check_mpi_configuration(num_threads)

    # Make sure the output either exists with --overwrite or --restart or
    # does not exist
    if os.path.isdir(calc_dir):
        if overwrite:
            rmtree(calc_dir)

        if (not restart) and (not overwrite):
            err = f"'{previous_run_dir}/{calc_dir}' already exists. Either remove\n"
            err += "it, specify --overwrite, or specify --restart.\n\n"
            os.chdir("..")
            raise FileExistsError(err)

    if os.path.isdir(calc_dir) and restart:

        if not os.path.isdir(os.path.join(calc_dir,"working","replicates")):
            err = "\ncould not restart the calculation. Please delete the directory\n"
            err += f"'{previous_run_dir}/{calc_dir}' and try again.\n\n"
            os.chdir("..")
            raise ValueError(err)

        # Load existing supervisor for the directory we are restarting. 
        supervisor = Supervisor(calc_dir)
        supervisor.update("calc_status","running")
        supervisor.event("Restarting calculation.")

        reconcile_bootstrap(supervisor.df,
                            supervisor.model,
                            supervisor.gene_tree,
                            supervisor.species_tree,
                            supervisor.reconciled_tree,
                            supervisor.run_parameters["allow_horizontal_transfer"],
                            supervisor.seed,
                            bootstrap_directory=None,
                            converge_cutoff=converge_cutoff,
                            restart="replicates",
                            overwrite=False,
                            supervisor=supervisor,
                            num_threads=num_threads,
                            threads_per_rep=threads_per_replicate,
                            generax_binary=generax_binary,
                            raxml_binary=raxml_binary)

    else:

        topiary.reconcile(prev_calculation=bootstrap_directory,
                          calc_dir=calc_dir,
                          bootstrap=True,
                          converge_cutoff=converge_cutoff,
                          overwrite=False,
                          num_threads=num_threads,
                          threads_per_rep=threads_per_replicate,
                          raxml_binary=raxml_binary,
                          generax_binary=generax_binary)

    os.chdir('..')

    # Create an html report for the calculation
    pipeline_report(pipeline_directory=previous_run_dir,
                    output_directory=os.path.join(previous_run_dir,"results"),
                    overwrite=True)