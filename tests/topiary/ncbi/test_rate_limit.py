"""
Tests for the NCBI request throttle.

NCBI returns HTTP 429 if queried faster than 3 requests/second (10 with an API
key). Every Entrez call site has to go through this, or a run with many species
gets throttled -- which is exactly what happened in CI.
"""

import time

import pytest

import topiary.ncbi


@pytest.fixture(autouse=True)
def reset_throttle():
    """
    Clear the module-level timestamp so each test starts from a clean state and
    does not inherit a sleep debt from its neighbour.
    """

    topiary.ncbi._last_request_time = 0.0
    yield
    topiary.ncbi._last_request_time = 0.0


def test_rate_limit_waits_between_back_to_back_calls(mocker):

    slept = []
    mocker.patch("topiary.ncbi.time.sleep", side_effect=slept.append)

    # Two calls with no real time in between: the second has to wait for
    # (almost) the whole interval.
    topiary.ncbi.rate_limit()
    topiary.ncbi.rate_limit()

    assert len(slept) == 1
    assert slept[0] == pytest.approx(topiary.ncbi.NCBI_REQUEST_FREQ, rel=0.1)


def test_rate_limit_does_not_wait_when_the_caller_was_already_slow(mocker):

    slept = []
    mocker.patch("topiary.ncbi.time.sleep", side_effect=slept.append)

    # Pretend the last request was long ago. No wait is owed.
    topiary.ncbi._last_request_time = time.monotonic() - 10.0

    topiary.ncbi.rate_limit()

    assert slept == []


def test_rate_limit_waits_only_for_the_remaining_interval(mocker):

    slept = []
    mocker.patch("topiary.ncbi.time.sleep", side_effect=slept.append)

    # Half the interval has already elapsed, so only half is owed. A bare
    # time.sleep(NCBI_REQUEST_FREQ) would sleep twice as long as needed.
    half = topiary.ncbi.NCBI_REQUEST_FREQ / 2
    topiary.ncbi._last_request_time = time.monotonic() - half

    topiary.ncbi.rate_limit()

    assert len(slept) == 1
    assert slept[0] == pytest.approx(half, abs=0.05)
    assert slept[0] < topiary.ncbi.NCBI_REQUEST_FREQ


def test_rate_limit_records_the_time_of_the_request(mocker):

    mocker.patch("topiary.ncbi.time.sleep")

    before = time.monotonic()
    topiary.ncbi.rate_limit()

    assert topiary.ncbi._last_request_time >= before


def test_rate_limit_actually_spaces_real_calls():
    """
    No mocking: three back-to-back calls must take at least two intervals of
    wall-clock time.
    """

    start = time.monotonic()
    for _ in range(3):
        topiary.ncbi.rate_limit()
    elapsed = time.monotonic() - start

    # First call is free; the next two each wait one interval.
    assert elapsed >= 2 * topiary.ncbi.NCBI_REQUEST_FREQ * 0.9


# ---------------------------------------------------------------------------
# Every Entrez call site must go through the throttle. These are the ones that
# were missing it and produced HTTP 429 in CI.
# ---------------------------------------------------------------------------

def test_get_taxid_is_throttled(mocker):

    from topiary.ncbi.entrez.taxid import get_taxid

    limiter = mocker.patch("topiary.ncbi.rate_limit")
    handle = mocker.MagicMock()
    mocker.patch("topiary.ncbi.entrez.taxid.Entrez.esearch", return_value=handle)
    mocker.patch("topiary.ncbi.entrez.taxid.Entrez.read",
                 return_value={"Count": "1", "IdList": ["9606"]})

    get_taxid("Homo sapiens")

    limiter.assert_called()


def test_get_mrca_taxid_is_throttled(mocker):

    from topiary.ncbi.entrez.mrca import get_mrca_taxid

    limiter = mocker.patch("topiary.ncbi.rate_limit")
    mocker.patch("topiary.ncbi.get_taxid", return_value=["9606", "10090"])
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.efetch")
    mocker.patch("topiary.ncbi.entrez.mrca.Entrez.read",
                 return_value=[{"LineageEx": [{"TaxId": "7711"}], "TaxId": "9606"},
                               {"LineageEx": [{"TaxId": "7711"}], "TaxId": "10090"}])

    get_mrca_taxid(["Homo sapiens", "Mus musculus"])

    limiter.assert_called()


def test__get_records_is_throttled(mocker):

    from topiary.ncbi.entrez.proteome import _get_records

    limiter = mocker.patch("topiary.ncbi.rate_limit")
    mocker.patch("topiary.ncbi.entrez.proteome.Entrez.esummary")
    mocker.patch("topiary.ncbi.entrez.proteome.Entrez.read",
                 return_value={"DocumentSummarySet":
                               {"DocumentSummary": [{"SpeciesName": "Homo sapiens"}]}})

    _get_records(["11968211"])

    limiter.assert_called()


def test_get_proteome_ids_is_throttled_once_per_query(mocker):
    """
    The esearch call sits inside a loop over filters, so it can fire more than
    once per call -- each pass has to be throttled, not just the first.
    """

    from topiary.ncbi.entrez.proteome import get_proteome_ids

    limiter = mocker.patch("topiary.ncbi.rate_limit")
    # proteome.py resolves the species through topiary.ncbi.get_taxid at call
    # time rather than importing it, so patch it there.
    mocker.patch("topiary.ncbi.get_taxid", return_value=9606)
    mocker.patch("topiary.ncbi.entrez.proteome.Entrez.esearch")

    # No ids come back, so the loop runs through every filter
    mocker.patch("topiary.ncbi.entrez.proteome.Entrez.read",
                 return_value={"IdList": []})

    get_proteome_ids(species="Homo sapiens")

    assert limiter.call_count >= 2
