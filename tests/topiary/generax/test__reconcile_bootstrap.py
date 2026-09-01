import pytest
import topiary

from topiary.generax._reconcile_bootstrap import _build_replicate_dirs
from topiary.generax._reconcile_bootstrap import _get_timeout_config
from topiary.generax._reconcile_bootstrap import _compute_replicate_timeout
from topiary.generax._reconcile_bootstrap import _terminate_process
from topiary.generax._reconcile_bootstrap import _launch_replicate
from topiary.generax._reconcile_bootstrap import _DEFAULT_TIMEOUT_CONFIG
from topiary.generax._generax import GENERAX_BINARY

import ete4 as ete

import os
import sys
import glob
import copy
import time
import subprocess


def test__launch_replicate_env(tmpdir, monkeypatch):

    monkeypatch.chdir(tmpdir)

    # The environment handed to the subprocess is the current process
    # environment (topiary no longer rewrites it for MPI).
    monkeypatch.setenv("TOPIARY_MARKER", "present")

    cmd = [sys.executable, "-c",
           "import os; print(os.environ.get('TOPIARY_MARKER', 'missing'))"]
    returncode, timed_out = _launch_replicate(cmd, timeout=30,
                                              stdout_path="stdout.log",
                                              stderr_path="stderr.log")
    assert returncode == 0
    assert timed_out is False
    with open("stdout.log") as f:
        assert "present" in f.read()


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
                "floor":1.0}
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


def test__terminate_process(tmpdir,monkeypatch):

    monkeypatch.chdir(tmpdir)

    # A long-running process is terminated. We launch it in the SAME session as
    # the test (no start_new_session) -- terminating it must not require a
    # process-group kill.
    proc = subprocess.Popen([sys.executable,"-c","import time; time.sleep(60)"])
    time.sleep(0.2)
    assert proc.poll() is None

    _terminate_process(proc)
    proc.wait(timeout=10)

    # The process is dead
    assert proc.poll() is not None


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

