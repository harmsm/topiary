import pytest

import topiary

import os
import shutil

@pytest.mark.run_generax
@pytest.mark.run_raxml
def test_integrated_minimal_ali_to_anc(tiny_phylo,tmpdir):
    """
    Full test of ali_to_anc pipeline without any complexity of arg checking
    etc. Goal is to catch major problems or changes to core functionality
    linking pieces together. We pass in a known species tree to
    avoid any interaction with open tree of life.
    """


    # Hangs dramatically on github workflows. Works fine locally, probably
    # because of challenge of configuring MPI remotely. (Not awesome) solution
    # is to mark as slow and not run slow tests on gh.

    df = tiny_phylo["initial-input/dataframe.csv"]
    species_tree = tiny_phylo["initial-input/species-tree.newick"]

    def _check_out_files(tiny_phylo,out_dir):

        key = "/".join([out_dir,"output","*"])
        expected_files = tiny_phylo[key]
        for e in expected_files:
            assert os.path.exists(os.path.join(out_dir,"output",
                                               os.path.basename(e)))


    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # find best model

    topiary.find_best_model(df=df,
                            calc_dir="00_find-best-model",
                            seed=12345,
                            model_matrices=["LG","JTT"],
                            model_rates=None,
                            model_freqs=None,
                            model_invariant=None)

    _check_out_files(tiny_phylo,"00_find-best-model")

    # -------------------------------------------------------------------------
    # generate ml tree

    topiary.generate_ml_tree(prev_calculation="00_find-best-model",
                             calc_dir="01_gene-tree")

    _check_out_files(tiny_phylo,"01_gene-tree")

    # -------------------------------------------------------------------------
    # reconcile (single-core generax by default)
    topiary.reconcile(prev_calculation="01_gene-tree",
                      calc_dir="02_reconcile",
                      species_tree=species_tree)

    _check_out_files(tiny_phylo,"02_reconcile")

    # -------------------------------------------------------------------------
    # infer ancestors

    topiary.generate_ancestors(prev_calculation="02_reconcile",
                               calc_dir="03_ancestors")

    _check_out_files(tiny_phylo,"03_ancestors")

    # -------------------------------------------------------------------------
    # gene tree bootstraps

    topiary.generate_bootstraps(prev_calculation="03_ancestors",
                                calc_dir="04_bootstraps",
                                num_threads=-1)

    _check_out_files(tiny_phylo,"04_bootstraps")

    # -------------------------------------------------------------------------
    # reconcile bootstraps (the crawler pipeline). We use the canned toy
    # bootstrap output as input rather than the tiny real bootstraps above, so
    # run it in its own directory. A single crawler invocation sets up, runs
    # every replicate, and aggregates.

    bs_run = "bootstrap-reconcile-run"
    os.mkdir(bs_run)
    shutil.copytree(tiny_phylo["04_bootstraps_toy"],
                    os.path.join(bs_run, "04_gene-tree-bootstraps"))

    topiary.bootstrap_reconcile(bs_run)

    # Because we are using the toy bootstraps rather than real bootstraps from
    # the last step, only check for the primary expected output.
    assert os.path.isfile(os.path.join(bs_run,
                                       "05_reconciled-tree-bootstraps",
                                       "output",
                                       "reconciled-tree_supports.newick"))

