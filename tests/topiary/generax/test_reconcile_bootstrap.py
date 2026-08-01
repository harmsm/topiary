import pytest
import topiary.generax._reconcile_bootstrap as rb
import os

def test__construct_args(mocker):
    # Mock mpi.get_hosts
    mock_get_hosts = mocker.patch("topiary.generax._reconcile_bootstrap.mpi.get_hosts")
    
    # Case 1: Simple case, factor of node size
    mock_get_hosts.return_value = ["n1", "n1", "n1", "n1"]
    kwargs_list, num_threads = rb._construct_args("dir", 0.03, 4, 2)
    
    assert num_threads == 2
    assert len(kwargs_list) == 2
    assert kwargs_list[0]["hosts"] == ["n1", "n1"]
    assert kwargs_list[1]["hosts"] == ["n1", "n1"]
    assert kwargs_list[0]["is_manager"] == True
    assert kwargs_list[1]["is_manager"] == False
    
    # Case 2: Spanning nodes (should print warning)
    mock_get_hosts.return_value = ["n1", "n1", "n2", "n2"]
    # We expect a warning because [n1, n1] is fine, but if threads_per_rep was 3:
    # Slice 1: [n1, n1, n2] -> warning
    
    # Capture stdout to check for warning
    import io
    from contextlib import redirect_stdout
    f = io.StringIO()
    with redirect_stdout(f):
        kwargs_list, num_threads = rb._construct_args("dir", 0.03, 4, 3)
    
    output = f.getvalue()
    assert "WARNING: A bootstrap replicate calculation spans multiple compute" in output
    assert "nodes (n1, n2)" in output or "nodes (n2, n1)" in output
    
    assert len(kwargs_list) == 2 # 4 slots / 3 per rep = 2 groups (one with 3, one with 1)
    assert kwargs_list[0]["hosts"] == ["n1", "n1", "n2"]
    assert kwargs_list[1]["hosts"] == ["n2"]

def test__generax_thread_function(mocker, tmpdir):
    # The subprocess launch now happens inside _launch_replicate; mock that seam
    # and verify the mpirun command is constructed with the right hosts.

    rep_dir = os.path.join(tmpdir, "replicates")
    os.makedirs(os.path.join(rep_dir, "00001"))
    with open(os.path.join(rep_dir, "00001", "run_generax.sh"), "w") as f:
        f.write("generax --args &> topiary.log\n")

    mocker.patch("topiary.generax._reconcile_bootstrap.mpi._get_mpi_oversubscribe",
                 return_value=False)
    mocker.patch("topiary.generax._reconcile_bootstrap.mpi.get_mpi_flags",
                 return_value=[])

    captured = {}
    def fake_launch(cmd, timeout, stdout_path, stderr_path):
        captured["cmd"] = cmd
        # Create a result tree so the replicate "succeeds".
        os.makedirs(os.path.join("result", "results", "reconcile"), exist_ok=True)
        with open(os.path.join("result", "results", "reconcile", "geneTree.newick"), "w") as f:
            f.write("((a,b),c);\n")
        return 0, False

    mocker.patch("topiary.generax._reconcile_bootstrap._launch_replicate",
                 side_effect=fake_launch)

    original_dir = os.getcwd()
    try:
        rb._generax_thread_function(rep_dir, 0.03, False, ["n1", "n1"], lock=None)
    finally:
        os.chdir(original_dir)

    # mpirun command constructed with the correct hosts
    cmd = captured["cmd"]
    assert "mpirun" in cmd
    assert "--host" in cmd
    assert "n1,n1" in cmd

    # Replicate ran to completion
    assert os.path.isfile(os.path.join(rep_dir, "00001", "completed"))
    assert not os.path.exists(os.path.join(rep_dir, "00001", "running"))

def test_generax_thread_function_failure(mocker, tmpdir):
    # An error out of the launch must be caught and turned into a `failed`
    # replicate (never left `running`).
    mocker.patch("topiary.generax._reconcile_bootstrap._launch_replicate",
                 side_effect=RuntimeError("MPI Crash!"))

    # Create a mock directory structure in tmpdir
    rep_dir = os.path.join(tmpdir, "replicates")
    os.makedirs(rep_dir)
    os.makedirs(os.path.join(rep_dir, "00001"))
    with open(os.path.join(rep_dir, "00001", "run_generax.sh"), "w") as f:
        f.write("generax --args &> topiary.log\n")

    # Mock mpi._get_mpi_oversubscribe and get_mpi_flags
    mocker.patch("topiary.generax._reconcile_bootstrap.mpi._get_mpi_oversubscribe",
                 return_value=False)
    mocker.patch("topiary.generax._reconcile_bootstrap.mpi.get_mpi_flags",
                 return_value=[])

    # Run thread function
    original_dir = os.getcwd()
    try:
        rb._generax_thread_function(rep_dir, 0.03, False, ["localhost"], lock=None)
    finally:
        os.chdir(original_dir)

    # Check if 'failed' file was created
    assert os.path.isfile(os.path.join(rep_dir, "00001", "failed"))
    # Check if 'running' file was removed (or never existed/was cleaned up)
    assert not os.path.exists(os.path.join(rep_dir, "00001", "running"))
