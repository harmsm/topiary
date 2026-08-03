import pytest
import topiary

import topiary.generax._reconcile_bootstrap as _rb
from topiary.generax._reconcile_bootstrap import _progress_bar
from topiary.generax._reconcile_bootstrap import _check_convergence
from topiary.generax._reconcile_bootstrap import _generax_thread_function
from topiary.generax._reconcile_bootstrap import _build_replicate_dirs
from topiary.generax._reconcile_bootstrap import _clean_replicate_dir
from topiary.generax._reconcile_bootstrap import _construct_args
from topiary.generax._reconcile_bootstrap import _run_bootstrap_calculations
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary.generax._reconcile_bootstrap import _LocalValue
from topiary.generax._reconcile_bootstrap import _get_timeout_config
from topiary.generax._reconcile_bootstrap import _compute_replicate_timeout
from topiary.generax._reconcile_bootstrap import _should_abort
from topiary.generax._reconcile_bootstrap import _kill_process_group
from topiary.generax._reconcile_bootstrap import _launch_replicate
from topiary.generax._reconcile_bootstrap import _DEFAULT_TIMEOUT_CONFIG
from topiary.generax._generax import GENERAX_BINARY
from topiary.raxml import RAXML_BINARY
from topiary._private import Supervisor
from topiary._private import mpi
from topiary._private.threads import MockLock

import ete4 as ete

import pandas as pd

import os
import sys
import glob
import shutil
import copy
import pathlib
import signal
import subprocess
import time
import multiprocessing as mp

def test__progress_bar(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Test just running...

    shutil.copytree(template,"test0")

    status_bar = mp.Process(target=_progress_bar,args=("test0",))
    status_bar.start()

    time.sleep(2)

    # Should not have completed.
    try:
        assert status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    # -------------------------------------------------------------------------
    # Test convergence -- suddenly write "skipped" to a bunch of directories

    shutil.copytree(template,"test1")

    status_bar = mp.Process(target=_progress_bar,args=("test1",))
    status_bar.start()

    for g in glob.glob(os.path.join("test1","0*")):

        # Turn running into completed
        if os.path.isfile(os.path.join(g,"running")):
            os.remove(os.path.join(g,"running"))
            pathlib.Path(os.path.join(g,"completed")).touch()

        # Put skipped in all directories that are not complete
        if not os.path.isfile(os.path.join(g,"completed")):
            pathlib.Path(os.path.join(g,"skipped")).touch()

    # Add timeout loop in case it takes a moment to finish
    start_time = time.time()
    while (time.time() - start_time) < 6:
        if not status_bar.is_alive():
            break
        time.sleep(0.2)

    # Should have completed because we wrote "skipped" and "complete" into all
    # directories
    try:
        assert not status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    # -------------------------------------------------------------------------
    # Test convergence -- write "completed" to all directories

    shutil.copytree(template,"test2")

    status_bar = mp.Process(target=_progress_bar,args=("test2",))
    status_bar.start()

    # Put completed in all directories
    for g in glob.glob(os.path.join("test2","0*")):

        # Clean up running so it's strictly completed
        if os.path.isfile(os.path.join(g,"running")):
            os.remove(os.path.join(g,"running"))

        # Turn running into completed
        if not os.path.isfile(os.path.join(g,"completed")):
            pathlib.Path(os.path.join(g,"completed")).touch()

    # Add timeout loop in case it takes a moment to finish
    status_bar.join(timeout=6)

    # Should have completed because we wrote "skipped" and "complete" into all
    # directories
    try:
        assert not status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    os.chdir(current_dir)

@pytest.mark.run_raxml
def test__check_convergence(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Converged

    test_dir = "test0"
    shutil.copytree(template,test_dir)
    c, df = _check_convergence(test_dir,converge_cutoff=0.5)
    assert c
    assert issubclass(type(df),pd.DataFrame)
    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # -------------------------------------------------------------------------
    # Not converged

    test_dir = "test1"
    shutil.copytree(template,test_dir)
    shutil.move(os.path.join(test_dir,"bs-trees_not-converged.newick"),
                os.path.join(test_dir,"bs-trees.newick"))

    c, df = _check_convergence(test_dir,converge_cutoff=0)
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    os.chdir(current_dir)

@pytest.mark.run_raxml
@pytest.mark.run_generax
def test__generax_thread_function(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"converge_cutoff":0.5,
                       "is_manager":False,
                       "lock":None}

    # --------------------------------------------------------------------------
    # Should run on last three, skipping first that are completed and running.

    test_dir = "test0"
    shutil.copytree(template,test_dir)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    out = _generax_thread_function(**kwargs)
    assert out is None

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))


    # --------------------------------------------------------------------------
    # Run test again, this time sending in two hosts. Behavior should be
    # identical.

    # There is a problem with generax where it screws up running in parallel on 
    # a fast processor with the small test dataset. (Generax expects a file will
    # be completely written out on one thread, but it's not quite finished
    # writing out when the other thread grabs it.) Loop through the test runner
    # 10 times to make sure it works at least once. Not perfect, but...

    def _test_runner():

        test_dir = "test1"
        if os.path.isdir(test_dir):
            shutil.rmtree(test_dir)

        shutil.copytree(template,test_dir)

        kwargs = copy.deepcopy(kwargs_template)
        kwargs["replicate_dir"] = test_dir
        kwargs["hosts"] = mpi.get_hosts(2)
        out = _generax_thread_function(**kwargs)
        assert out is None

        for i in range(10):
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        for i in range(10,12):
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        for i in range(12,15):
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        # Make sure calculation did *not* run in 00001, which should not have a
        # results directory because it had a 'completed' file.
        assert not os.path.exists(os.path.join(test_dir,"00001","results"))
        assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))

    for i in range(10):
        try:
            _test_runner()
            break
        except RuntimeError:
            continue

    # --------------------------------------------------------------------------
    # Run test again with a single host and manager = True. Should run once, then
    # make last two directories have "skipped" because calculation has converged

    test_dir = "test2"
    shutil.copytree(template,test_dir)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    kwargs["is_manager"] = True
    out = _generax_thread_function(**kwargs)
    assert out[0]
    assert issubclass(type(out[1]),pd.DataFrame)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,13):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(13,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))


    # --------------------------------------------------------------------------
    # Run test gain with a single host and manager = True, but unconverged
    # bootstrap replicates. Should complete calculations and return converged
    # False.

    test_dir = "test3"
    shutil.copytree(template,test_dir)
    shutil.move(os.path.join(test_dir,"bs-trees_not-converged.newick"),
                os.path.join(test_dir,"bs-trees.newick"))

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    kwargs["is_manager"] = True
    kwargs["converge_cutoff"] = 0.0000001
    out = _generax_thread_function(**kwargs)
    assert not out[0]
    assert issubclass(type(out[1]),pd.DataFrame)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))

    os.chdir(current_dir)

@pytest.mark.run_generax
def test__build_replicate_dirs(small_phylo,tmpdir):

    input_dir = small_phylo["05_gene-tree-bootstraps_toy/output/"]
    df = topiary.read_dataframe(small_phylo["05_gene-tree-bootstraps_toy/input/dataframe.csv"])

    f = open(small_phylo["model.txt"])
    model = f.read().strip()
    f.close()

    gene_tree = small_phylo["final-output/gene-tree.newick"]
    species_tree = small_phylo["final-output/species-tree.newick"]
    bootstrap_directory = small_phylo["05_gene-tree-bootstraps_toy/output/bootstrap_replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"df":df,
                       "model":model,
                       "gene_tree":gene_tree,
                       "species_tree":"species_tree",
                       "allow_horizontal_transfer":True,
                       "seed":12345,
                       "bootstrap_directory":bootstrap_directory,
                       "overwrite":False,
                       "generax_binary":GENERAX_BINARY}

    os.mkdir("test0")
    os.chdir("test0")
    kwargs = copy.deepcopy(kwargs_template)

    _build_replicate_dirs(**kwargs)

    # Make sure that directories are made with correct files
    expected_files = set(["alignment.phy",
                          "control.txt",
                          "gene_tree.newick",
                          "mapping.link",
                          "run_generax.sh",
                          "species_tree.newick"])
    expected = "replicates"
    assert len(list(glob.glob(os.path.join(expected,"0*")))) == 4
    for g in glob.glob(os.path.join(expected,"0*")):
        dirs = set(os.listdir(g))
        assert dirs == expected_files

    # Make sure that the correct files are identical between references and
    # that the correct files are different
    def _read_file(some_file):
        out = []
        with open(some_file,'r') as f:
            for line in f:
                out.append(line.strip())
        return "\n".join(out)

    # Read in alignments
    alignments = []
    for i in range(4):
        alignments.append(_read_file(os.path.join(input_dir,"bootstrap_replicates",f"bsmsa_000{i+1}.phy")))

    trees = []
    with open(os.path.join(input_dir,"bootstrap_replicates","bs-trees.newick")) as f:
        for line in f:
            T = ete.Tree(line.strip(),parser=0)
            trees.append(T)

    ref_check = {}
    should_be_same = ["control.txt","mapping.link","run_generax.sh","species_tree.newick"]

    dir_list = glob.glob(os.path.join(expected,"0*"))
    dir_list.sort()
    for i, g in enumerate(dir_list):

        this_check = {}
        for f in glob.glob(os.path.join(g,"*")):
            this_check[os.path.basename(f)] = _read_file(f)

            if i == 0:
                ref_check[os.path.basename(f)] = _read_file(f)

        assert this_check["alignment.phy"] == alignments[i]

        # Make sure the gene trees are correctly brought in
        T = ete.Tree(this_check["gene_tree.newick"],parser=0)
        assert T.robinson_foulds(trees[i],unrooted_trees=True)[0] == 0

        for k in this_check:
            if k in should_be_same:
                assert ref_check[k] == this_check[k]

    os.chdir(current_dir)

def test__clean_replicate_dir(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    test_dir = "test0"
    shutil.copytree(template,test_dir)
    pathlib.Path(os.path.join(test_dir,"00013","skipped")).touch()

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,13):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(13,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    _clean_replicate_dir(test_dir)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))


    os.chdir(current_dir)

@pytest.mark.run_generax
def test__construct_args():

    kwargs_template = {"replicate_dir":"test",
                       "converge_cutoff":0.03,
                       "num_threads":1,
                       "threads_per_rep":1}

    # 1 thread, 1 per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 1
    assert num_threads == 1
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.03
    assert len(kwargs_list[0]["hosts"]) == 1
    assert kwargs_list[0]["is_manager"]

    # 2 threads, 1 per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 2
    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 2
    assert num_threads == 2
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.03
    assert len(kwargs_list[0]["hosts"]) == 1
    assert kwargs_list[0]["is_manager"]

    assert kwargs_list[1]["replicate_dir"] == "test"
    assert kwargs_list[1]["converge_cutoff"] == 0.03
    assert len(kwargs_list[1]["hosts"]) == 1
    assert not kwargs_list[1]["is_manager"]

    # 2 threads, 2 threads per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 2
    kwargs["threads_per_rep"] = 2
    kwargs["converge_cutoff"] = 0.5

    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 1
    assert num_threads == 1
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.5
    assert len(kwargs_list[0]["hosts"]) == 2
    assert kwargs_list[0]["is_manager"]


@pytest.mark.run_generax
def test__run_bootstrap_calculations(small_phylo,tmpdir):

    # This basically makes sure that a calculation is run in every directory.
    # It does not check quality/correctness of output

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    shutil.copytree(small_phylo["06_reconciled-tree-bootstraps_toy/working/replicates"],
                    "replicates_template")

    kwargs_template = {"replicate_dir":"replicates",
                       "converge_cutoff":0.03,
                       "num_threads":1,
                       "threads_per_rep":1}

    # Single thread
    rep_dir = "test0"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))

    # Run multithreaded

    rep_dir = "test1"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir
    kwargs["num_threads"] = 2

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))


    # Run multithreaded, multiple slots

    rep_dir = "test2"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir
    kwargs["num_threads"] = 2
    kwargs["threads_per_rep"] = 2

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))

    os.chdir(current_dir)


@pytest.mark.run_raxml
@pytest.mark.run_generax
def test_reconcile_bootstrap(small_phylo,tmpdir):

    df_csv = small_phylo["initial-input/dataframe.csv"]
    df = topiary.read_dataframe(df_csv)
    gene_tree = small_phylo["final-output/gene-tree.newick"]
    species_tree = small_phylo["initial-input/species-tree.newick"]
    reconciled_tree = small_phylo["final-output/reconciled-tree.newick"]
    input_bootstrap_directory = small_phylo["05_gene-tree-bootstraps_toy/output/bootstrap_replicates/"]
    f = open(small_phylo["model.txt"],"r")
    model = f.read().strip()
    f.close()

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"df":df,
                       "model":model,
                       "gene_tree":gene_tree,
                       "species_tree":species_tree,
                       "reconciled_tree":reconciled_tree,
                       "allow_horizontal_transfer":True,
                       "bootstrap_directory":input_bootstrap_directory,
                       "converge_cutoff":0.03,
                       "seed":True,
                       "restart":None,
                       "overwrite":False,
                       "supervisor":None,
                       "num_threads":1,
                       "threads_per_rep":1,
                       "generax_binary":GENERAX_BINARY,
                       "raxml_binary":RAXML_BINARY}

    supervisor = Supervisor()
    supervisor.create_calc_dir("test0",
                               calc_type="test0",
                               df=df,
                               gene_tree=gene_tree,
                               model=model)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["supervisor"] = supervisor

    tT = reconcile_bootstrap(**kwargs)

    output_dir = supervisor.output_dir
    expected_files = ["dataframe.csv",
                      "reconciliations.txt",
                      "summary-tree.pdf",
                      "reconciled-tree_supports.newick"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(output_dir,f))

    new_T = ete.Tree(os.path.join(output_dir,"reconciled-tree_supports.newick"),parser=0)
    old_T = ete.Tree(reconciled_tree)

    # Topology should *not* have changed
    assert new_T.robinson_foulds(old_T,unrooted_trees=True)[0] == 0

    # Make sure it now has supports
    for n in new_T.traverse():
        if not n.is_leaf:
            print(n.support)

    os.chdir(current_dir)


# -----------------------------------------------------------------------------
# Tests for the per-replicate timeout / failure-handling machinery.
# -----------------------------------------------------------------------------

def test_mocklock_acquire():

    # MockLock.acquire should mirror the multiprocessing lock proxy interface:
    # accept blocking/timeout and return True.
    lock = MockLock()
    assert lock.acquire() is True
    assert lock.acquire(timeout=5) is True
    assert lock.acquire(blocking=False) is True
    assert lock.acquire(True,10) is True
    assert lock.release() is None


def test__LocalValue():

    v = _LocalValue()
    assert v.value == 0

    v = _LocalValue(7)
    assert v.value == 7

    v.value += 3
    assert v.value == 10


def test__get_timeout_config():

    # None -> a copy of the defaults (not the same object)
    config = _get_timeout_config(None)
    assert config == _DEFAULT_TIMEOUT_CONFIG
    assert config is not _DEFAULT_TIMEOUT_CONFIG

    # Partial override merges on top of defaults
    config = _get_timeout_config({"factor":10.0})
    assert config["factor"] == 10.0
    assert config["ceiling"] == _DEFAULT_TIMEOUT_CONFIG["ceiling"]
    assert config["floor"] == _DEFAULT_TIMEOUT_CONFIG["floor"]

    # Full override
    override = {"factor":2.0,
                "ceiling":10.0,
                "floor":1.0,
                "max_failed_fraction":0.5,
                "max_failed_floor":1}
    config = _get_timeout_config(override)
    assert config == override

    # Unrecognized key raises
    with pytest.raises(ValueError):
        _get_timeout_config({"not_a_key":1})


def test__compute_replicate_timeout():

    config = {"factor":3.0,"ceiling":1000.0,"floor":10.0,
              "max_failed_fraction":0.1,"max_failed_floor":5}

    # No sample threshold -> always the ceiling, even with data
    assert _compute_replicate_timeout([],None,config) == 1000.0
    assert _compute_replicate_timeout([1,2,3],None,config) == 1000.0

    # Not enough samples yet -> ceiling
    assert _compute_replicate_timeout([],3,config) == 1000.0
    assert _compute_replicate_timeout([5.0,5.0],3,config) == 1000.0

    # Enough samples, factor dominates the floor
    #   max = 20, factor*max = 60 > floor (10) -> 60
    assert _compute_replicate_timeout([10.0,20.0,5.0],3,config) == 60.0

    # Enough samples, floor dominates
    #   max = 1, factor*max = 3 < floor (10) -> 10
    assert _compute_replicate_timeout([1.0,1.0,1.0],3,config) == 10.0


def test__should_abort():

    config = {"factor":3.0,"ceiling":1000.0,"floor":10.0,
              "max_failed_fraction":0.1,"max_failed_floor":5}

    # Below the failure floor -> never abort, regardless of fraction
    assert _should_abort(0,100,config) is False
    assert _should_abort(4,10,config) is False

    # num_total is None or non-positive -> never abort
    assert _should_abort(100,None,config) is False
    assert _should_abort(100,0,config) is False

    # At/above floor and above the fraction -> abort
    #   6 > 0.1 * 10 = 1.0 -> True
    assert _should_abort(6,10,config) is True

    # At/above floor but not above the fraction -> do not abort
    #   6 failures out of 100 -> 6 > 10.0 is False
    assert _should_abort(6,100,config) is False

    # Exactly equal to the fraction is not "more than" -> do not abort
    #   10 == 0.1 * 100 -> not > -> False
    assert _should_abort(10,100,config) is False


def test__kill_process_group(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    # Launch a shell that spawns a sleeping grandchild, all in a new session so
    # they form their own process group.
    proc = subprocess.Popen(["sh","-c","sleep 60"],start_new_session=True)

    # Let the group come up and grab its id.
    time.sleep(0.3)
    pgid = os.getpgid(proc.pid)

    # signal 0 -> group currently exists
    os.killpg(pgid,0)

    _kill_process_group(proc)
    proc.wait(timeout=10)

    # The direct child is dead
    assert proc.poll() is not None

    # The whole group (including the sleep grandchild) is gone
    time.sleep(0.3)
    with pytest.raises(ProcessLookupError):
        os.killpg(pgid,0)


def test__launch_replicate(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Normal completion: stdout/stderr captured to files, returncode reported,
    # not timed out.

    cmd = [sys.executable,"-c",
           "import sys; sys.stdout.write('hello out'); sys.stderr.write('hello err')"]
    returncode, timed_out = _launch_replicate(cmd,
                                              timeout=30,
                                              stdout_path="stdout.log",
                                              stderr_path="stderr.log")
    assert returncode == 0
    assert timed_out is False

    with open("stdout.log") as f:
        assert "hello out" in f.read()
    with open("stderr.log") as f:
        assert "hello err" in f.read()

    # -------------------------------------------------------------------------
    # Non-zero exit code is reported (but not treated as a timeout)

    cmd = [sys.executable,"-c","import sys; sys.exit(3)"]
    returncode, timed_out = _launch_replicate(cmd,
                                              timeout=30,
                                              stdout_path="stdout2.log",
                                              stderr_path="stderr2.log")
    assert returncode == 3
    assert timed_out is False

    # -------------------------------------------------------------------------
    # Timeout: a long-running process is killed and flagged as timed out. The
    # call must return well before the process would have finished on its own.

    cmd = [sys.executable,"-c","import time; time.sleep(60)"]
    start = time.time()
    returncode, timed_out = _launch_replicate(cmd,
                                              timeout=0.5,
                                              stdout_path="stdout3.log",
                                              stderr_path="stderr3.log")
    elapsed = time.time() - start
    assert timed_out is True
    assert elapsed < 30


def test__launch_replicate_env(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    # The environment handed to the subprocess must come from mpi.get_mpi_env.
    monkeypatch.setattr(_rb.mpi,"get_mpi_env",
                        lambda: {"TOPIARY_MARKER":"present"})

    cmd = [sys.executable,"-c",
           "import os; print(os.environ.get('TOPIARY_MARKER','missing'))"]
    returncode, timed_out = _launch_replicate(cmd,
                                              timeout=30,
                                              stdout_path="stdout.log",
                                              stderr_path="stderr.log")
    assert returncode == 0
    assert timed_out is False
    with open("stdout.log") as f:
        assert "present" in f.read()


def _build_fake_replicate_dir(repdir,names):
    """
    Build a minimal replicate directory tree with `run_generax.sh` files for
    each replicate in `names`.
    """

    os.mkdir(repdir)
    for name in names:
        d = os.path.join(repdir,name)
        os.mkdir(d)
        with open(os.path.join(d,"run_generax.sh"),"w") as f:
            f.write("generax --families control.txt &> topiary.log\n")


def _make_fake_launch(plan):
    """
    Build a fake `_launch_replicate` that behaves according to `plan`, a dict
    mapping replicate directory name -> one of "success", "timeout", "notree".
    """

    result_tree = os.path.join("result","results","reconcile","geneTree.newick")

    def fake_launch(cmd,timeout,stdout_path,stderr_path):

        # We are inside the replicate directory when this is called.
        with open(stdout_path,"w") as f:
            f.write("fake stdout\n")
        with open(stderr_path,"w") as f:
            f.write("fake stderr\n")

        outcome = plan[os.path.basename(os.getcwd())]

        if outcome == "success":
            os.makedirs(os.path.dirname(result_tree),exist_ok=True)
            with open(result_tree,"w") as f:
                f.write("(A:1,B:1);\n")
            return 0, False

        if outcome == "timeout":
            return None, True

        if outcome == "raise":
            raise RuntimeError("boom")

        # "notree" -- exited (non-zero) without producing a result tree
        return 1, False

    return fake_launch


def test__generax_thread_function_failure_handling(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    names = ["00001","00002","00003","00004","00005"]
    _build_fake_replicate_dir("replicates",names)

    # 00002 times out, 00003 produces no tree; the rest succeed.
    plan = {"00001":"success",
            "00002":"timeout",
            "00003":"notree",
            "00004":"success",
            "00005":"success"}
    monkeypatch.setattr(_rb,"_launch_replicate",_make_fake_launch(plan))

    durations = []
    fail_count = _LocalValue(0)

    out = _generax_thread_function(replicate_dir="replicates",
                                   converge_cutoff=0.5,
                                   is_manager=False,
                                   hosts=["localhost"],
                                   lock=None,
                                   durations=durations,
                                   fail_count=fail_count,
                                   total_replicates=len(names),
                                   sample_threshold=None,
                                   timeout_config=None)
    assert out is None

    # Successful replicates are marked completed (not failed)
    for name in ["00001","00004","00005"]:
        assert os.path.isfile(os.path.join("replicates",name,"completed"))
        assert not os.path.isfile(os.path.join("replicates",name,"failed"))
        assert not os.path.isfile(os.path.join("replicates",name,"running"))

    # Failed replicates (timeout + missing tree) are marked failed (not completed)
    for name in ["00002","00003"]:
        assert os.path.isfile(os.path.join("replicates",name,"failed"))
        assert not os.path.isfile(os.path.join("replicates",name,"completed"))
        assert not os.path.isfile(os.path.join("replicates",name,"running"))

    # Runtime recorded only for the three successes
    assert len(durations) == 3

    # Failure counter incremented for the two failures
    assert fail_count.value == 2

    # bs-trees.newick has exactly the three successful trees
    with open(os.path.join("replicates","bs-trees.newick")) as f:
        lines = [line for line in f if line.strip() != ""]
    assert len(lines) == 3

    # Failure reason annotated into the stderr log
    with open(os.path.join("replicates","00002","stderr.log")) as f:
        assert "timeout" in f.read()
    with open(os.path.join("replicates","00003","stderr.log")) as f:
        assert "result tree" in f.read()


def test__generax_thread_function_unexpected_exception(tmpdir,monkeypatch):

    # An unexpected exception inside the launch must be caught and turned into a
    # normal replicate failure (directory flagged `failed`, not left `running`),
    # so it can never wedge the whole calculation.

    monkeypatch.chdir(tmpdir)

    names = ["00001","00002"]
    _build_fake_replicate_dir("replicates",names)

    plan = {"00001":"raise","00002":"success"}
    monkeypatch.setattr(_rb,"_launch_replicate",_make_fake_launch(plan))

    durations = []
    fail_count = _LocalValue(0)

    out = _generax_thread_function(replicate_dir="replicates",
                                   converge_cutoff=0.5,
                                   is_manager=False,
                                   hosts=["localhost"],
                                   lock=None,
                                   durations=durations,
                                   fail_count=fail_count,
                                   total_replicates=len(names),
                                   sample_threshold=None,
                                   timeout_config=None)
    assert out is None

    # The exception replicate is failed (and not stuck running)
    assert os.path.isfile(os.path.join("replicates","00001","failed"))
    assert not os.path.isfile(os.path.join("replicates","00001","running"))
    assert not os.path.isfile(os.path.join("replicates","00001","completed"))

    # The following replicate still ran to completion
    assert os.path.isfile(os.path.join("replicates","00002","completed"))

    assert fail_count.value == 1
    with open(os.path.join("replicates","00001","stderr.log")) as f:
        assert "unexpected error" in f.read()


def test__generax_thread_function_circuit_breaker(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    names = ["00001","00002","00003","00004","00005"]
    _build_fake_replicate_dir("replicates",names)
    repdir = os.path.abspath("replicates")

    # Everything fails.
    plan = {name:"notree" for name in names}
    monkeypatch.setattr(_rb,"_launch_replicate",_make_fake_launch(plan))

    durations = []
    fail_count = _LocalValue(0)

    # floor of 2 failures, 10% fraction of 5 replicates -> abort on the 2nd
    # failure (2 >= 2 and 2 > 0.5).
    timeout_config = {"max_failed_floor":2,"max_failed_fraction":0.1}

    with pytest.raises(RuntimeError):
        _generax_thread_function(replicate_dir="replicates",
                                 converge_cutoff=0.5,
                                 is_manager=False,
                                 hosts=["localhost"],
                                 lock=None,
                                 durations=durations,
                                 fail_count=fail_count,
                                 total_replicates=len(names),
                                 sample_threshold=None,
                                 timeout_config=timeout_config)

    # Aborting mid-run can leave the working directory changed; restore it so
    # relative-path bookkeeping in later tests is unaffected.
    monkeypatch.chdir(tmpdir)

    # Aborted after the second failure; later replicates never ran.
    assert fail_count.value == 2
    assert os.path.isfile(os.path.join(repdir,"00001","failed"))
    assert os.path.isfile(os.path.join(repdir,"00002","failed"))
    assert not os.path.isfile(os.path.join(repdir,"00003","failed"))


def test__generax_thread_function_adaptive_timeout(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    names = ["00001","00002","00003"]
    _build_fake_replicate_dir("replicates",names)

    plan = {name:"success" for name in names}

    # Record the timeout handed to each replicate so we can verify the adaptive
    # behavior. The first replicate should get the (long) ceiling; once we have
    # a sample, subsequent replicates get factor*max(durations), floored.
    seen_timeouts = []
    base_fake = _make_fake_launch(plan)

    def recording_fake(cmd,timeout,stdout_path,stderr_path):
        seen_timeouts.append(timeout)
        return base_fake(cmd,timeout,stdout_path,stderr_path)

    monkeypatch.setattr(_rb,"_launch_replicate",recording_fake)

    durations = []
    fail_count = _LocalValue(0)
    timeout_config = {"factor":3.0,"ceiling":9999.0,"floor":1.0,
                      "max_failed_fraction":0.1,"max_failed_floor":5}

    _generax_thread_function(replicate_dir="replicates",
                             converge_cutoff=0.5,
                             is_manager=False,
                             hosts=["localhost"],
                             lock=None,
                             durations=durations,
                             fail_count=fail_count,
                             total_replicates=len(names),
                             sample_threshold=1,
                             timeout_config=timeout_config)

    # First replicate: no samples yet -> ceiling.
    assert seen_timeouts[0] == 9999.0

    # Subsequent replicates: adaptive. With a sample_threshold of 1, once one
    # replicate has completed we switch to factor*max(durations), floored at 1.
    assert len(seen_timeouts) == 3
    for t in seen_timeouts[1:]:
        assert t != 9999.0
        assert t >= 1.0
