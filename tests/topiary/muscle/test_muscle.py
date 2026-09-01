import pytest

import os

import pandas as pd

import topiary
from topiary.muscle.muscle import _run_muscle
from topiary.muscle.muscle import align


class _FakePopen:
    """
    Stand-in for subprocess.Popen as _run_muscle uses it: iterate stdout, close
    it, then wait() for a return code.
    """

    def __init__(self, lines=("aligning 50.0%", "aligning 100.0%"), return_code=0):
        self.stdout = iter(lines)
        self._return_code = return_code
        self.closed = False
        self.waited = False

    def close_stdout(self):
        self.closed = True

    def wait(self):
        self.waited = True
        return self._return_code


def _install_fake_popen(mocker, popen, capture):
    """
    Patch subprocess.Popen inside the muscle module, recording the command.
    """

    def _factory(cmd, **kwargs):
        capture.append(cmd)
        popen.stdout = iter(list(popen.stdout))
        # Popen.stdout.close() is called directly by _run_muscle
        popen.stdout = _ClosableIter(popen.stdout, popen)
        return popen

    return mocker.patch("topiary.muscle.muscle.subprocess.Popen", side_effect=_factory)


class _ClosableIter:
    def __init__(self, it, owner):
        self._it = it
        self._owner = owner

    def __iter__(self):
        return self._it

    def close(self):
        self._owner.closed = True


def _muscle_version(mocker, version):
    return mocker.patch("topiary.muscle.muscle.installed.check_muscle",
                        return_value=("muscle", version, None))


def test__run_muscle_builds_muscle5_command(mocker):

    _muscle_version(mocker, (5, 1, 0))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(), cmds)

    _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")

    assert cmds[0] == ["muscle", "-align", "in.fasta", "-output", "out.fasta"]


def test__run_muscle_super5_swaps_the_algorithm(mocker):

    _muscle_version(mocker, (5, 1, 0))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(), cmds)

    _run_muscle("in.fasta", "out.fasta", True, True, [], "muscle")

    assert cmds[0][1] == "-super5"


def test__run_muscle_builds_legacy_command_for_muscle3(mocker):

    _muscle_version(mocker, (3, 8, 1551))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(), cmds)

    _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")

    assert cmds[0] == ["muscle", "-in", "in.fasta", "-out", "out.fasta"]


def test__run_muscle_appends_user_arguments(mocker):

    _muscle_version(mocker, (5, 1, 0))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(), cmds)

    _run_muscle("in.fasta", "out.fasta", False, True,
                ["-replicates", "20"], "muscle")

    assert cmds[0][-2:] == ["-replicates", "20"]


def test__run_muscle_raises_when_muscle_is_missing(mocker):

    _muscle_version(mocker, (-2, -2, -2))

    with pytest.raises(RuntimeError):
        _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")


def test__run_muscle_raises_when_muscle_is_crashing(mocker):

    mocker.patch("topiary.muscle.muscle.installed.check_muscle",
                 return_value=("muscle", (-1, -1, -1), {"stderr": "boom"}))

    with pytest.raises(RuntimeError):
        _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")


def test__run_muscle_warns_but_proceeds_on_unknown_version(mocker):

    _muscle_version(mocker, (0, 0, 0))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(), cmds)

    with pytest.warns(UserWarning):
        _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")

    # Unknown version is assumed to be >= 5, so we get the modern command line
    assert cmds[0][1] == "-align"


def test__run_muscle_raises_on_nonzero_return(mocker):

    _muscle_version(mocker, (5, 1, 0))
    cmds = []
    _install_fake_popen(mocker, _FakePopen(return_code=1), cmds)

    import subprocess
    with pytest.raises(subprocess.CalledProcessError):
        _run_muscle("in.fasta", "out.fasta", False, True, [], "muscle")


def test_align_rejects_bad_input_type():

    with pytest.raises(ValueError):
        align(1.0)

    with pytest.raises(ValueError):
        align(["not", "a", "dataframe"])


def test_align_fasta_input_requires_output_fasta(tmpdir):

    fasta = os.path.join(tmpdir, "in.fasta")
    with open(fasta, "w") as f:
        f.write(">seq1\nMASTPD\n")

    # Exists, but no output specified
    with pytest.raises(ValueError):
        align(fasta, output_fasta=None)


def test_align_raises_on_missing_fasta():

    with pytest.raises(FileNotFoundError):
        align("this-file-does-not-exist.fasta", output_fasta="out.fasta")


def test_align_fasta_input_calls_muscle(mocker, tmpdir):

    fasta = os.path.join(tmpdir, "in.fasta")
    with open(fasta, "w") as f:
        f.write(">seq1\nMASTPD\n")

    run = mocker.patch("topiary.muscle.muscle._run_muscle")

    out = align(fasta, output_fasta="out.fasta")

    # Aligning a file returns None; the result is the file muscle wrote
    assert out is None
    run.assert_called_once()
    assert run.call_args[0][0] == fasta
    assert run.call_args[0][1] == "out.fasta"


def test_align_dataframe_returns_dataframe_and_cleans_up(mocker, tmpdir):

    os.chdir(tmpdir)

    df = pd.DataFrame({"name": ["A", "B"],
                       "species": ["Homo sapiens", "Mus musculus"],
                       "sequence": ["MAST", "MASS"],
                       "keep": [True, True],
                       "uid": ["aaaaaaaaaa", "bbbbbbbbbb"]})

    def _fake_run_muscle(input_fasta, output_fasta, *args, **kwargs):
        # muscle would write an aligned fasta; do the same so read_fasta_into
        # has something to read.
        with open(output_fasta, "w") as f:
            f.write(">aaaaaaaaaa\nMAST\n>bbbbbbbbbb\nMASS\n")

    mocker.patch("topiary.muscle.muscle._run_muscle",
                 side_effect=_fake_run_muscle)

    out = align(df, silent=True)

    assert isinstance(out, pd.DataFrame)
    assert "alignment" in out.columns

    # Both the temporary input and the temporary output are cleaned up
    leftover = [f for f in os.listdir(".") if f.startswith("topiary-tmp_")]
    assert leftover == []


def test_align_dataframe_keeps_named_output_file(mocker, tmpdir):

    os.chdir(tmpdir)

    df = pd.DataFrame({"name": ["A"],
                       "species": ["Homo sapiens"],
                       "sequence": ["MAST"],
                       "keep": [True],
                       "uid": ["aaaaaaaaaa"]})

    def _fake_run_muscle(input_fasta, output_fasta, *args, **kwargs):
        with open(output_fasta, "w") as f:
            f.write(">aaaaaaaaaa\nMAST\n")

    mocker.patch("topiary.muscle.muscle._run_muscle",
                 side_effect=_fake_run_muscle)

    align(df, output_fasta="kept.fasta", silent=True)

    # A user-specified output file is not treated as temporary
    assert os.path.isfile("kept.fasta")


def test_align_cleans_up_temp_input_when_muscle_fails(mocker, tmpdir):
    """
    A muscle failure must not leave topiary-tmp_* scratch files behind in the
    caller's working directory.
    """

    os.chdir(tmpdir)

    df = pd.DataFrame({"name": ["A"],
                       "species": ["Homo sapiens"],
                       "sequence": ["MAST"],
                       "keep": [True],
                       "uid": ["aaaaaaaaaa"]})

    mocker.patch("topiary.muscle.muscle._run_muscle",
                 side_effect=RuntimeError("muscle fell over"))

    with pytest.raises(RuntimeError):
        align(df, silent=True)

    leftover = [f for f in os.listdir(".") if f.startswith("topiary-tmp_")]
    assert leftover == []
