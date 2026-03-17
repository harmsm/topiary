import pytest
import topiary.generax._reconcile_bootstrap as rb
from unittest.mock import MagicMock

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

def test__generax_thread_function(mocker):
    # Mock subprocess.run
    mock_run = mocker.patch("topiary.generax._reconcile_bootstrap.subprocess.run")
    mock_run.return_value = MagicMock(returncode=0, stdout=b"out", stderr=b"err")
    
    # Mock other things to avoid file IO
    mocker.patch("os.chdir")
    mocker.patch("os.listdir", return_value=["00001"])
    mocker.patch("os.path.isdir", return_value=True)
    mocker.patch("os.path.isfile", side_effect=[False, False, False, True, True]) # completed, running, skipped, run_generax.sh, result_tree
    mocker.patch("pathlib.Path.touch")
    mocker.patch("os.remove")
    
    # Define a custom side-effect for open to handle different files
    def side_effect_open(filename, mode='r'):
        mock = MagicMock()
        if 'r' in mode:
            if "run_generax.sh" in filename:
                mock.__enter__.return_value.read.return_value = "generax --args\n"
            elif "geneTree.newick" in filename:
                mock.__enter__.return_value.read.return_value = "((a,b),c);"
            else:
                mock.__enter__.return_value.read.return_value = ""
        return mock

    mocker.patch("builtins.open", side_effect=side_effect_open)
    
    # Mock mpi._get_mpi_oversubscribe
    mocker.patch("topiary.generax._reconcile_bootstrap.mpi._get_mpi_oversubscribe", return_value=False)
    mocker.patch("topiary.generax._reconcile_bootstrap.get_mpi_env", return_value={})

    # Mock lock
    mock_lock = MagicMock()
    
    # Run thread function
    rb._generax_thread_function("replicate_dir", 0.03, False, ["n1", "n1"], lock=mock_lock)
    
    # Check if mpirun was called with correct hosts
    found_mpirun = False
    for call_args in mock_run.call_args_list:
        args = call_args[0][0]
        if "mpirun" in args:
            found_mpirun = True
            assert "--host" in args
            assert "n1,n1" in args
    
    assert found_mpirun
