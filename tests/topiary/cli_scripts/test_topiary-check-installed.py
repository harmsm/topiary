import pytest


import sys
import os
import subprocess

from topiary.cli_scripts.check_installed import main

@pytest.mark.smoke
def test_main():

    # simple test. Make sure it runs. We test validation stack directly
    # elsewhere.

    # Get location of binary
    location = os.path.dirname(os.path.realpath(__file__))
    test_bin = "topiary-check-installed"

    if os.name == "nt":
        base_cmd = [sys.executable,test_bin]
    else:
        base_cmd = [test_bin]

    # The console script has to exist on the PATH and exit cleanly.
    ret = subprocess.run(base_cmd,capture_output=True)

    assert ret.returncode == 0

    # And it has to actually report on the stack rather than printing nothing
    assert len(ret.stdout) > 0


def test_main_surfaces_crash(mocker,capsys):
    """
    When validate_stack reports a crashing binary (e.g. an architecture-
    incompatible raxml-ng dying with 'Illegal instruction'), main() should
    print that detailed message rather than silently assuming a $PATH problem.
    """

    import topiary

    crash_message = ("raxml-ng (/some/path/raxml-ng): killed by signal SIGILL\n"
                     "Illegal instruction (core dumped)\n"
                     "This usually means the binary is incompatible with this "
                     "machine (architecture).")

    mocker.patch.object(topiary._private.installed,
                        "validate_stack",
                        side_effect=RuntimeError(crash_message))

    # Should not raise -- the RuntimeError is caught and reported.
    main()

    captured = capsys.readouterr()
    assert "SIGILL" in captured.out
    assert "Illegal instruction (core dumped)" in captured.out
    assert "$PATH" in captured.out
