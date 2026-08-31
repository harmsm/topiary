
import pytest

from topiary._private.installed import _version_checker
from topiary._private.installed import _build_diagnostic
from topiary._private.installed import _format_diagnostic
from topiary._private.installed import check_git
from topiary._private.installed import check_muscle
from topiary._private.installed import check_generax
from topiary._private.installed import check_raxml
from topiary._private.installed import check_blastp
from topiary._private.installed import check_makeblastdb
from topiary._private.installed import _compare_versions
from topiary._private.installed import validate_stack

from topiary.generax import GENERAX_BINARY

import warnings
import os
import sys
import signal

class _FakeReturn:
    """Minimal stand-in for subprocess.CompletedProcess."""
    def __init__(self,returncode,stdout=b"",stderr=b""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr

def test__version_checker():

    # Good check -- should work
    cmd = ["git","--version"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v, d = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert len(v) > 1
    assert int(v[0]) > 0
    assert d is None

    # Should fail with binary not found (-2,-2,-2)
    cmd = ["not_really_a_binary"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v, d = _version_checker(cmd,_version_slicer)
    assert b is None
    assert type(v) is tuple
    assert v == (-2,-2,-2)
    assert d is None

    # Should fail with could not run (-1,-1,-1). A binary that exits non-zero
    # should also come back with a populated diagnostic (returncode != 0, no
    # signal, captured output).
    cmd = ["git","--bad_argument_not_recognized"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v, d = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert v == (-1,-1,-1)
    assert d is not None
    assert d["returncode"] != 0
    assert d["signal"] is None
    # something explaining the failure should have been captured
    assert (d["stderr"].strip() != "") or (d["stdout"].strip() != "")

    # Should fail with ran but could not get version (0,0,0)
    cmd = ["git","--version"]
    def _version_slicer(ret):
        # bad parsing call -- last split()[1] will throw an IndexError
        return ret.stdout.decode().split()[2].strip().split()[1]

    b, v, d = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert v == (0,0,0)
    assert d is None


@pytest.mark.skipif(os.name == "nt",
                    reason="signal-based process death is not portable to Windows")
def test__version_checker_signal():
    """
    A binary killed by a signal (e.g. the SIGILL / 'Illegal instruction' seen
    when running an architecture-incompatible binary) should be reported as
    (-1,-1,-1) with a diagnostic naming the signal.
    """

    # Small program that kills itself with SIGILL, mimicking an
    # architecture-incompatible binary crashing with "Illegal instruction".
    cmd = [sys.executable,"-c",
           "import os,signal; os.kill(os.getpid(), signal.SIGILL)"]

    def _version_slicer(ret):
        return ret.stdout.decode().strip()

    b, v, d = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert v == (-1,-1,-1)
    assert d is not None
    assert d["returncode"] == -int(signal.SIGILL)
    assert d["signal"] == "SIGILL"


def test__build_diagnostic():

    # Non-zero exit, no signal
    d = _build_diagnostic(_FakeReturn(1,stdout=b"out",stderr=b"boom"))
    assert d["returncode"] == 1
    assert d["signal"] is None
    assert d["stdout"] == "out"
    assert d["stderr"] == "boom"

    # Killed by a signal -- returncode is negative and signal name resolved
    d = _build_diagnostic(_FakeReturn(-int(signal.SIGILL)))
    assert d["returncode"] == -int(signal.SIGILL)
    assert d["signal"] == "SIGILL"

    # Non-decodable / None streams should not blow up
    d = _build_diagnostic(_FakeReturn(1,stdout=None,stderr=None))
    assert d["stdout"] == ""
    assert d["stderr"] == ""


def test__format_diagnostic():

    assert _format_diagnostic(None) == "found but did not run"

    d = {"returncode":1,"signal":None,"stdout":"","stderr":""}
    assert "return code" in _format_diagnostic(d)
    assert "1" in _format_diagnostic(d)

    d = {"returncode":-4,"signal":"SIGILL","stdout":"","stderr":""}
    assert "SIGILL" in _format_diagnostic(d)

def test_check_muscle():

    binary, version, diagnostic = check_muscle()

    if version == (-2,-2,-2):
        warnings.warn("muscle not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("muscle is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("muscle is installed but we cannot parse its version string!")

@pytest.mark.run_generax
def test_check_generax():

    binary, version, diagnostic = check_generax()

    if version == (-2,-2,-2):
        warnings.warn("generax not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("generax is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("generax is installed but we cannot parse its version string!")

@pytest.mark.run_raxml
def test_check_raxml():

    binary, version, diagnostic = check_raxml()

    if version == (-2,-2,-2):
        warnings.warn("raxml-ng not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("raxml-ng is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("raxml-ng is installed but we cannot parse its version string!")

def test_check_blastp():

    binary, version, diagnostic = check_blastp()

    if version == (-2,-2,-2):
        warnings.warn("blastp not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("blastp is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("blastp is installed but we cannot parse its version string!")

def test_check_makeblastdb():

    binary, version, diagnostic = check_makeblastdb()

    if version == (-2,-2,-2):
        warnings.warn("makeblastdb not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("makeblastdb is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("makeblastdb is installed but we cannot parse its version string!")

def test_check_git():

    binary, version, diagnostic = check_git()

    if version == (-2,-2,-2):
        warnings.warn("git not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("git is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("git is installed but we cannot parse its version string!")

# mark for run_generax because that's why we worry about MPI in the first place.
def test__compare_versions():

    # good version string; only specified element (1,) matches
    out = _compare_versions(("1","0"),(1,))
    assert out is True

    # bad version string; only specified element (1,) matches
    out = _compare_versions(("1","1b"),(1,))
    assert out is True

    # can't compare last element -- ambiguous
    out = _compare_versions(("1","1b"),(1,1))
    assert out is None

    # first comparable element bad -- should not pass
    out = _compare_versions(("1","1b"),(2,1))
    assert out is False

    # good version string; second element too low
    out = _compare_versions(("1","0"),(1,2))
    assert out is False

    # bad version string; can't check second position
    out = _compare_versions(("1","1b"),(1,2))
    assert out is None

    # First position too low
    out = _compare_versions(("0","1b"),(1,2))
    assert out is False

    # Shoudl pass -- last element high enough
    out = _compare_versions(("1","1","1"),(1,1,0))
    assert out is True

    # Should pass, matches version exactly
    out = _compare_versions(("1","1","1"),(1,1,1))
    assert out is True


def test_validate_stack():

    # not an amazing test, but at least checks core logic of whether or not
    # version is high enough.

    validate_stack([{"program":"git",
                               "min_version":(0,0,1),
                               "must_pass":True}])

    with pytest.raises(RuntimeError):
        validate_stack([{"program":"git",
                                   "min_version":(10000000,0,1),
                                   "must_pass":True}])


def test_validate_stack_crash(mocker):
    """
    A program that is present but crashes when run should raise a RuntimeError
    whose message names the crash (and signal) and points at an
    architecture/compilation problem -- not merely a $PATH problem.
    """

    import topiary._private.installed as installed

    diagnostic = {"returncode":-int(signal.SIGILL),
                  "signal":"SIGILL",
                  "stdout":"",
                  "stderr":"Illegal instruction (core dumped)"}

    # Pretend raxml-ng is found at a path but crashed with SIGILL
    mocker.patch.object(installed,
                        "check_raxml",
                        return_value=("/some/path/raxml-ng",(-1,-1,-1),diagnostic))

    with pytest.raises(RuntimeError) as exc_info:
        validate_stack([{"program":"raxml-ng",
                         "min_version":(1,1),
                         "must_pass":True}])

    msg = str(exc_info.value)
    assert "raxml-ng" in msg
    assert "SIGILL" in msg
    assert "Illegal instruction (core dumped)" in msg
    # should point the user at the architecture/compilation cause
    assert "architecture" in msg.lower()


def test_validate_stack_not_found(mocker):
    """
    A program that is simply missing from the $PATH should raise a RuntimeError
    that talks about the $PATH, and should NOT be reported as a crash.
    """

    import topiary._private.installed as installed

    mocker.patch.object(installed,
                        "check_raxml",
                        return_value=(None,(-2,-2,-2),None))

    with pytest.raises(RuntimeError) as exc_info:
        validate_stack([{"program":"raxml-ng",
                         "min_version":(1,1),
                         "must_pass":True}])

    msg = str(exc_info.value)
    assert "raxml-ng" in msg
    assert "$PATH" in msg
    # not-found is not a crash, so the crash-only guidance should be absent
    assert "crashed when" not in msg
