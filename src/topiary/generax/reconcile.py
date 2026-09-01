"""
Reconcile a gene tree with a species tree using generax.
"""

import topiary
from topiary._private import Supervisor
from topiary._private import check

from ._reconcile_no_bootstrap import reconcile_no_bootstrap

from ._generax import GENERAX_BINARY
from topiary.raxml import RAXML_BINARY

import subprocess
import os

def reconcile(prev_calculation=None,
              df=None,
              model=None,
              gene_tree=None,
              species_tree=None,
              reconciled_tree=None,
              allow_horizontal_transfer=None,
              seed=None,
              bootstrap=False,
              calc_dir="reconcile",
              overwrite=False,
              generax_launch="",
              generax_binary=GENERAX_BINARY,
              raxml_binary=RAXML_BINARY):
    """
    Reconcile the gene tree to the species tree using generax.

    Parameters
    ----------
    prev_calculation : str or Supervisor, optional
        previously completed calculation. Should either be a directory
        containing the calculation (e.g. the directory with run_parameters.json,
        input, working, output) or a Supervisor instance with a calculation
        loaded. Function will load dataframe, model, gene_tree, and
        reconciled_tree from the previous run. If this is not specified, `df`,
        `model`, `gene_tree` and `reconciled_tree` arguments must be
        specified.
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
        names; all other data will be dropped. NOTE: this is required if
        bootstrap = True.
    allow_horizontal_transfer : bool, optional
        whether to allow horizontal transfer during reconciliation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model. If None, use
        whatever is in prev_calculation. If this is not specified, default to
        True.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap : bool, default=False
        deprecated. Bootstrap reconciliation is now run via the
        topiary-bootstrap-reconcile crawler; passing True here raises an error.
    calc_dir: str, default="reconcile"
        name of calc_dir directory
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    supervisor : Supervisor, optional
        supervisor instance to keep track of inputs and outputs
    generax_launch : str, default=""
        launcher prefix prepended to the generax command (e.g. "mpirun -np 8").
        Empty string runs generax as a single process. The caller owns the
        launcher and the resources it needs.
    generax_binary : str, optional
        what generax binary to use
    raxml_binary : str, optional
        what raxml binary to use

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

    # Make sure that raxml is in the path
    try:
        subprocess.run([raxml_binary],capture_output=True)
    except FileNotFoundError:
        err = f"\nraxml binary '{raxml_binary}' not found in path\n\n"
        raise ValueError(err)

    # Bootstrap reconciliation moved to its own re-entrant, MPI-free crawler.
    if bootstrap:
        err = "\nbootstrap reconciliation is no longer run through reconcile().\n"
        err += "Use the topiary-bootstrap-reconcile command (a filesystem\n"
        err += "crawler) instead.\n\n"
        raise ValueError(err)

    # --------------------------------------------------------------------------
    # Load/parse calculation inputs

    # Load in previous calculation. Three possibilities here: prev_calculation
    # is a supervisor (just use it); prev_calculation is a directory (create a
    # supervisor from it); prev_calculation is None (create an empty supervisor).
    if isinstance(prev_calculation,Supervisor):
        supervisor = prev_calculation
    else:
        supervisor = Supervisor(calc_dir=prev_calculation)

    # Create a calculation directory
    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type="reconcile",
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree,
                               species_tree=species_tree,
                               reconciled_tree=reconciled_tree,
                               model=model)

    if allow_horizontal_transfer is None:
        if "allow_horizontal_transfer" not in supervisor.run_parameters:
            supervisor.run_parameters["allow_horizontal_transfer"] = True
    else:
        allow_horizontal_transfer = check.check_bool(allow_horizontal_transfer)
        supervisor.update("allow_horizontal_transfer",allow_horizontal_transfer)

    supervisor.check_required(required_values=["model","allow_horizontal_transfer"],
                              required_files=["alignment.phy","dataframe.csv",
                                              "gene-tree.newick"])

    # If no species tree is given, infer and load into supervisor
    if supervisor.species_tree is None:
        species_tree, dropped = topiary.df_to_species_tree(supervisor.df)
        species_tree_out = os.path.join(supervisor.input_dir,"species-tree.newick")
        species_tree.write(outfile=species_tree_out,parser=5)
        supervisor.update("species_tree",species_tree_out)

    allow_ht = supervisor.run_parameters["allow_horizontal_transfer"]

    return reconcile_no_bootstrap(df=supervisor.df,
                                  model=supervisor.model,
                                  gene_tree=supervisor.gene_tree,
                                  species_tree=supervisor.species_tree,
                                  allow_horizontal_transfer=allow_ht,
                                  seed=supervisor.seed,
                                  overwrite=overwrite,
                                  supervisor=supervisor,
                                  generax_launch=generax_launch,
                                  generax_binary=generax_binary)
