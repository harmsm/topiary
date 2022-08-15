"""
Reconcile a gene tree with a species tree using generax.
"""

import topiary
from topiary._private import Supervisor
from topiary._private.mpi import get_num_slots
from topiary._private.mpi import check_mpi_configuration

from ._reconcile_bootstrap import reconcile_bootstrap
from ._reconcile_no_bootstrap import reconcile_no_bootstrap

from ._generax import GENERAX_BINARY

import subprocess
import os

def reconcile(previous_dir=None,
              df=None,
              model=None,
              gene_tree=None,
              species_tree=None,
              reconciled_tree=None,
              allow_horizontal_transfer=True,
              bootstrap=False,
              calc_dir="reconcile",
              overwrite=False,
              supervisor=None,
              num_threads=-1,
              generax_binary=GENERAX_BINARY):
    """
    Reconcile the gene tree to the species tree using generax.

    Parameters
    ----------
    previous_dir : str, optional
        directory containing previous calculation. function will grab the the
        csv, model, and tree from the previous run. If this is not specified,
        `df`, `model`, and `tree_file` arguments must be specified.
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `previous_dir` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `previous_dir`
        if specified.
    gene_tree : str, ete3.Tree, dendropy.tree, optional
        gene tree file for calculation. Will override tree in `previous_dir`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names and branch lengths; all other data will be dropped.
    species_tree : str, ete3.Tree, dendropy.tree, optional
        species tree file for calculation. Will override tree in `previous_dir`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped.
    reconciled_tree : str, ete3.Tree, dendropy.tree, optional
        reconciled tree file for calculation. Will override tree in `previous_dir`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped. NOTE: this is required if
        bootstrap = True.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    bootstrap: bool, default=False
        whether or not to do bootstrap replicates. if True, previous_dir must
        point to a raxml ml_bootstrap run
    calc_dir: str, default="reconcile"
        name of calc_dir directory
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    supervisor : Supervisor, optional
        supervisor instance to keep track of inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1 use all available.
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    # Make sure that generax is in the path
    try:
        subprocess.run([generax_binary],capture_output=True)
    except FileNotFoundError:
        err = f"\ngenerax binary '{generax_binary}' not found in path\n\n"
        raise ValueError(err)

    # Get number of slots
    if num_threads == -1:
        num_threads = get_num_slots(generax_binary)

    # Check sanity of mpi configuration/number of slots
    check_mpi_configuration(num_threads,generax_binary)

    # --------------------------------------------------------------------------
    # Load/parse calculation inputs

    # Load existing
    if supervisor is None:
        supervisor = Supervisor(previous_dir)

    # Create a calculation directory
    if bootstrap:
        calc_type = "reconcile_bootstrap"
    else:
        calc_type = "reconcile"

    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type=calc_type,
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree,
                               species_tree=species_tree,
                               reconciled_tree=reconciled_tree,
                               model=model)

    supervisor.update("allow_horizontal_transfer",allow_horizontal_transfer)

    supervisor.check_required(required_values=["model","allow_horizontal_transfer"],
                              required_files=["alignment.phy","dataframe.csv",
                                              "gene-tree.newick"])

    # If no species tree is given, infer and load into supervisor
    if supervisor.species_tree is None:
        species_tree, dropped = topiary.df_to_species_tree(supervisor.df)
        species_tree_out = os.path.join(supervisor.input_dir,"species-tree.newick")
        species_tree.write(outfile=species_tree_out,format=5)
        supervisor.update("species_tree",species_tree_out)

    allow_ht = supervisor.run_parameters["allow_horizontal_transfer"]

    if not bootstrap:

        return reconcile_no_bootstrap(df=supervisor.df,
                                      model=supervisor.model,
                                      gene_tree=supervisor.gene_tree,
                                      species_tree=supervisor.species_tree,
                                      allow_horizontal_transfer=allow_ht,
                                      overwrite=overwrite,
                                      supervisor=supervisor,
                                      num_threads=num_threads,
                                      generax_binary=generax_binary)

    else:

        try:
            prev_calc_type = supervisor.previous_entries[-1]["calc_type"]
        except:
            prev_calc_type = None

        if prev_calc_type != "ml_bootstrap":
            err = "\nboostrap reconciliation can only be started using a\n"
            err += "previous 'ml_bootstrap' run as its input.\n"
            raise ValueError(err)

        # Make sure bootstrap directory exists
        bs_dir = os.path.join(supervisor.previous_entries[-1]["calc_dir"],
                              "output","bootstrap_replicates")

        if not os.path.isdir(bs_dir):
            err = f"\ninput directory '{supervisor.previous_entries[-1]['calc_dir']}'\n"
            err += "does not have an output/bootstrap_replicates directory. Was\n"
            err += "this calculation run with bootstrap=True?\n\n"
            raise FileNotFoundError(err)

        # Make sure the supervisor has a reconciled tree loaded. This is needed
        # because this is what the bootstraps will be mapped onto. Generally
        # this will the ML reconciled tree from a previous reconciliation
        # calculation in the pipeline.
        supervisor.check_required(required_files=["reconciled-tree.newick"])

        return reconcile_bootstrap(df=supervisor.df,
                                   model=supervisor.model,
                                   gene_tree=supervisor.gene_tree,
                                   species_tree=supervisor.species_tree,
                                   reconciled_tree=supervisor.reconciled_tree,
                                   allow_horizontal_transfer=allow_ht,
                                   bootstrap_directory=bs_dir,
                                   overwrite=overwrite,
                                   supervisor=supervisor,
                                   num_threads=num_threads,
                                   generax_binary=generax_binary)