"""
Tests of the test harness itself -- the autouse fixtures in tests/conftest.py
that keep one test from corrupting the next.

These use pytest's `pytester` fixture, which runs a throwaway pytest session in
a temporary directory. That lets us assert on cross-test contamination without
making the tests in this file order-dependent.
"""

import os

import pytest


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
