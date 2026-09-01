"""
Regression test for the Open Tree of Life file-descriptor leak.

The opentree client records every API response in OT.ws.call_history, and each
record holds the requests.Response object, which keeps its socket open. Nothing
releases them, so a process leaked one descriptor per OTL call. Running the
whole test suite made several hundred such calls and exhausted the default
macOS limit of 256 descriptors, which is why the suite died partway through
with "Too many open files" while individual test files passed.

topiary disables the recording in topiary/opentree/util.py. These tests make
sure it stays disabled.
"""

import os

import pytest

from topiary.opentree.util import OT, species_to_ott


def _open_fd_count():
    for fd_dir in ("/proc/self/fd", "/dev/fd"):
        try:
            return len(os.listdir(fd_dir))
        except OSError:
            continue
    return None


def test_opentree_call_recording_is_disabled():
    """
    The cheap, hermetic version: just assert the setting importing topiary is
    supposed to have applied. Runs in the default suite, no network needed.
    """

    assert OT.ws._store_api_calls is False, (
        "OT.ws._store_api_calls is back on. Every Open Tree of Life call will "
        "pin its socket open and the suite will run out of file descriptors.")


@pytest.mark.network
def test_opentree_calls_do_not_leak_descriptors():
    """
    The real thing: hit the API repeatedly and confirm the descriptor count
    does not climb.
    """

    if _open_fd_count() is None:
        pytest.skip("cannot count file descriptors on this platform")

    # One warm-up call so any lazily-created connection/pool exists before we
    # take the baseline.
    species_to_ott(["Homo sapiens"])

    before = _open_fd_count()

    for _ in range(5):
        species_to_ott(["Homo sapiens"])

    after = _open_fd_count()

    assert after <= before, (
        f"leaked {after - before} file descriptors over 5 Open Tree of Life "
        f"calls ({before} -> {after})")
