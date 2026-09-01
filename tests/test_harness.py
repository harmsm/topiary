"""
Tests of the test harness itself -- the autouse fixtures in tests/conftest.py
that keep one test from corrupting the next.

These use pytest's `pytester` fixture, which runs a throwaway pytest session in
a temporary directory. That lets us assert on cross-test contamination without
making the tests in this file order-dependent.
"""

import os

import pytest


# The tier-assignment and subprocess-guard logic from tests/conftest.py, cut
# down to just the parts under test. Kept as a literal (rather than importing
# the real conftest) so pytester sessions are isolated from this session's
# options and markers.
_TIER_CONFTEST = '''
import pytest, subprocess

_OPT_IN_MARKERS = frozenset(["network","run_ncbi_server","run_blast",
                             "run_generax","run_raxml"])

def pytest_configure(config):
    for m in ("unit","smoke","integration","network","run_ncbi_server",
              "run_blast","run_generax","run_raxml"):
        config.addinivalue_line("markers", m + ": test marker")

def pytest_collection_modifyitems(config, items):
    for item in items:
        if len(set(item.keywords) & _OPT_IN_MARKERS) > 0:
            item.add_marker(pytest.mark.integration)
        elif "smoke" in item.keywords:
            pass
        else:
            item.add_marker(pytest.mark.unit)

class SubprocessBlocked(Exception):
    pass

@pytest.fixture(autouse=True)
def block_subprocess(request):
    if "unit" not in request.keywords:
        yield
        return
    real_popen = subprocess.Popen
    def guarded_popen(*args, **kwargs):
        cmd = args[0] if args else kwargs.get("args")
        raise SubprocessBlocked(
            f"{request.node.nodeid} is in the `unit` tier but tried to run {cmd!r}.")
    subprocess.Popen = guarded_popen
    try:
        yield
    finally:
        subprocess.Popen = real_popen
'''


def test_cwd_is_restored_after_a_test_that_fails_mid_chdir(pytester):
    """
    A test that chdirs away and then fails before its restore line must not
    leave the process in that directory.

    This is the failure mode that produced topiary's path-dependent test
    failures: the suite passed when tests ran individually and failed when they
    ran together, because one test's leaked cwd broke unrelated later tests.
    """

    pytester.makeconftest("""
        import pytest, os, tempfile

        @pytest.fixture(autouse=True)
        def restore_cwd():
            start_dir = os.getcwd()
            try:
                yield
            finally:
                try:
                    os.chdir(start_dir)
                except OSError:
                    os.chdir(tempfile.gettempdir())
    """)

    pytester.makepyfile("""
        import os

        def test_leaks_cwd(tmpdir):
            current_dir = os.getcwd()
            os.chdir(tmpdir)
            raise RuntimeError("fails before the restore")
            os.chdir(current_dir)

        def test_bystander(request):
            # Must still be in the directory pytest started from.
            assert os.getcwd() == str(request.config.invocation_params.dir)
    """)

    result = pytester.runpytest()

    # The genuinely-broken test still fails; the innocent one is unaffected.
    result.assert_outcomes(failed=1, passed=1)


def test_cwd_fixture_survives_its_start_directory_being_deleted(pytester):
    """
    If a test deletes the directory it started in, the fixture must still leave
    the process somewhere valid. Otherwise os.getcwd() raises for every
    subsequent test and the whole session collapses.
    """

    pytester.makeconftest("""
        import pytest, os, tempfile

        @pytest.fixture(autouse=True)
        def restore_cwd():
            start_dir = os.getcwd()
            try:
                yield
            finally:
                try:
                    os.chdir(start_dir)
                except OSError:
                    os.chdir(tempfile.gettempdir())
    """)

    pytester.makepyfile("""
        import os, shutil, tempfile

        def test_destroys_its_own_cwd():
            doomed = tempfile.mkdtemp()
            os.chdir(doomed)
            shutil.rmtree(doomed)

        def test_bystander():
            # Would raise FileNotFoundError if the fixture had left us inside
            # the deleted directory.
            assert os.path.isdir(os.getcwd())
    """)

    result = pytester.runpytest()
    result.assert_outcomes(passed=2)


def test_the_real_conftest_restores_cwd(tmp_path):
    """
    Belt and braces: the fixture above is a copy of the one in
    tests/conftest.py. Confirm the real one is actually installed and active in
    this session by chdir-ing away and letting teardown put us back.

    Paired with test_the_real_conftest_restored_cwd below, which checks the
    result. The pair is order-dependent by design -- pytest runs tests within a
    file in definition order -- and exists only to catch the case where someone
    deletes the autouse fixture but leaves the pytester tests above passing.
    """

    test_the_real_conftest_restores_cwd.start = os.getcwd()
    os.chdir(tmp_path)

    assert os.getcwd() != test_the_real_conftest_restores_cwd.start


def test_the_real_conftest_restored_cwd():

    expected = getattr(test_the_real_conftest_restores_cwd, "start", None)
    if expected is None:
        pytest.skip("companion test did not run")

    assert os.getcwd() == expected


def test_every_test_gets_exactly_one_tier(pytester):
    """
    Every collected test must land in exactly one of unit / smoke /
    integration. Two tiers on one test, or none, means the selection commands
    in CLAUDE.md silently run the wrong set.
    """

    pytester.makeconftest(_TIER_CONFTEST)

    pytester.makepyfile("""
        import pytest

        def test_plain():
            pass

        @pytest.mark.smoke
        def test_marked_smoke():
            pass

        @pytest.mark.run_generax
        def test_needs_generax():
            pass
    """)

    for tier in ("unit", "smoke", "integration"):
        result = pytester.runpytest("--collect-only", "-q", "-m", tier)
        result.stdout.fnmatch_lines([f"1/3 tests collected (2 deselected)*"])

    # And no test carries two tiers at once.
    for pair in ("unit and smoke", "unit and integration", "smoke and integration"):
        result = pytester.runpytest("--collect-only", "-q", "-m", pair)
        result.stdout.fnmatch_lines(["*no tests collected (3 deselected)*"])


def test_unit_tier_cannot_shell_out(pytester):
    """
    The `unit` tier claims to be hermetic. Enforce it: a unit test that spawns
    a subprocess must fail rather than quietly making the fast suite slow.
    """

    pytester.makeconftest(_TIER_CONFTEST)

    pytester.makepyfile("""
        import subprocess

        def test_sneaky_subprocess():
            subprocess.Popen(["echo", "hello"])
    """)

    result = pytester.runpytest()
    result.assert_outcomes(failed=1)
    result.stdout.fnmatch_lines(["*is in the `unit` tier but tried to run*"])


@pytest.mark.smoke
def test_smoke_tier_may_shell_out(pytester):
    """
    The counterpart: the guard must not fire outside the unit tier.

    Marked smoke itself: pytester runs the inner session in-process, so if this
    outer test were in the unit tier its own guard would patch subprocess.Popen
    and the inner test would fail for the wrong reason.
    """

    pytester.makeconftest(_TIER_CONFTEST)

    pytester.makepyfile("""
        import pytest, subprocess

        @pytest.mark.smoke
        def test_allowed_subprocess():
            p = subprocess.Popen(["echo", "hello"])
            p.wait()
    """)

    result = pytester.runpytest()
    result.assert_outcomes(passed=1)
