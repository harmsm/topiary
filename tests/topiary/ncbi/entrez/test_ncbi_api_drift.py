"""
Drift check against the live NCBI Entrez API.

Everything else in the suite mocks NCBI, which means a change to the shape of
NCBI's responses would not be noticed until a user's run broke. These tests make
a small number of real calls and assert on the *structure* topiary depends on --
the specific keys it reads out of each response -- rather than on values, which
change legitimately as NCBI updates its data.

Two rules keep this from becoming another flaky CI failure:

1. Every test asserts on shape, never on a value that NCBI is entitled to
   change (an assembly's update date, how many taxids match a name, and so on).
2. A transport-level problem -- HTTP 429, a 5xx, a timeout, a DNS failure -- is
   a *skip*, not a failure. NCBI being busy is not drift. Only a response that
   parses differently than topiary expects fails.

Run with `pytest tests --run-ncbi-server`. CI runs just this file, nightly, on a
single matrix configuration; see .github/workflows/python-app.yml.
"""

import datetime
import http.client
import re
import socket
import urllib.error

import pytest

import topiary
from topiary.ncbi.entrez.mrca import get_mrca_taxid
from topiary.ncbi.entrez.proteome import _get_records
from topiary.ncbi.entrez.proteome import get_proteome_ids
from topiary.ncbi.entrez.taxid import get_taxid


# Transport failures that mean "NCBI is unhappy right now", not "the API
# changed". These skip rather than fail.
TRANSIENT = (urllib.error.HTTPError,
             urllib.error.URLError,
             http.client.HTTPException,
             socket.timeout,
             TimeoutError,
             ConnectionError)


def _live(description, call):
    """
    Run a live NCBI call, turning a transport failure into a skip.

    Returns whatever `call` returns.
    """

    try:
        return call()
    except TRANSIENT as e:
        pytest.skip(f"NCBI unavailable for {description} ({type(e).__name__}: {e}); "
                    f"not a drift failure")


# A human assembly id, stable enough to look up. Used only to fetch a record --
# nothing is asserted about its contents beyond the fields topiary reads.
STABLE_ASSEMBLY_ID = "11968211"


@pytest.mark.run_ncbi_server
def test_taxonomy_esearch_still_has_the_fields_get_taxid_reads():
    """
    get_taxid reads record["Count"] and record["IdList"] off a taxonomy
    esearch, and treats Count as an int.
    """

    taxid = _live("taxonomy esearch", lambda: get_taxid("Homo sapiens"))

    # get_taxid returns the parsed taxid; if Count/IdList had moved or been
    # renamed this would have raised rather than returned.
    assert taxid is not None
    assert str(taxid).isdigit(), f"expected a numeric taxid, got {taxid!r}"


@pytest.mark.run_ncbi_server
def test_taxonomy_esearch_still_returns_errorlist_on_no_match():
    """
    When a name matches nothing, Count comes back different from the number of
    species asked about, and topiary reads record["ErrorList"] to build a
    useful message. If NCBI stopped returning ErrorList, topiary would raise a
    bare KeyError instead -- so the drift check is that the failure is a
    RuntimeError, not that it fails at all.

    Note "Not a species"-style phrases are not reliable here: NCBI's term
    parser matches some of them to the root taxon and returns a real id. Use a
    nonsense binomial.
    """

    with pytest.raises(RuntimeError):
        _live("taxonomy esearch (no match)",
              lambda: get_taxid(["Xyzzy plughwerp"]))


@pytest.mark.run_ncbi_server
def test_taxonomy_efetch_still_has_lineage_fields():
    """
    get_mrca_taxid walks record["LineageEx"], reading ["TaxId"] off each entry,
    then record["TaxId"] for the taxon itself.

    Assert on the response *shape*, not on the MRCA it computes. NCBI
    intermittently returns a record with an empty LineageEx under load, which
    makes get_mrca_taxid fall back to the root of life -- that is NCBI being
    busy, not a schema change, so it skips.
    """

    from Bio import Entrez

    taxids = _live("taxonomy esearch for efetch",
                   lambda: get_taxid(["Homo sapiens", "Mus musculus"]))

    def _fetch():
        topiary.ncbi.rate_limit()
        handle = Entrez.efetch(db="taxonomy", id=",".join(taxids), retmode="xml")
        try:
            return Entrez.read(handle)
        finally:
            handle.close()

    records = _live("taxonomy efetch", _fetch)

    assert len(records) == len(taxids)

    for record in records:

        # These two keys are the schema topiary depends on. Losing either is
        # drift and should fail.
        assert "TaxId" in record, (
            f"taxonomy efetch record no longer has 'TaxId'. "
            f"topiary reads it in ncbi/entrez/mrca.py. "
            f"Keys present: {sorted(record.keys())}")
        assert "LineageEx" in record, (
            f"taxonomy efetch record no longer has 'LineageEx'. "
            f"topiary reads it in ncbi/entrez/mrca.py. "
            f"Keys present: {sorted(record.keys())}")

        if len(record["LineageEx"]) == 0:
            pytest.skip("NCBI returned an empty LineageEx for a taxon that has "
                        "one; transient, not drift")

        for taxon in record["LineageEx"]:
            assert "TaxId" in taxon, (
                f"LineageEx entries no longer carry 'TaxId'. "
                f"Keys present: {sorted(taxon.keys())}")


@pytest.mark.run_ncbi_server
def test_get_mrca_taxid_returns_an_integer_taxid():
    """
    topiary's own contract on top of that response. Deliberately does not assert
    *which* ancestor -- that depends on NCBI's taxonomy content, which is
    allowed to change and is occasionally incomplete.
    """

    mrca = _live("get_mrca_taxid",
                 lambda: get_mrca_taxid(["Homo sapiens", "Mus musculus"]))

    assert isinstance(mrca, int)
    assert mrca >= 1


@pytest.mark.run_ncbi_server
def test_assembly_esearch_still_returns_an_idlist():
    """
    get_proteome_ids reads search_record["IdList"] off an assembly esearch.
    """

    ids, err = _live("assembly esearch",
                     lambda: get_proteome_ids(species="Homo sapiens"))

    assert err is None, f"assembly esearch failed: {err}"
    assert ids is not None
    assert len(ids) > 0
    assert all(str(i).strip() != "" for i in ids)


@pytest.mark.run_ncbi_server
def test_assembly_esummary_still_has_every_field_topiary_reads():
    """
    This is the richest drift surface. _get_records reaches into
    esummary_record["DocumentSummarySet"]["DocumentSummary"], and
    _get_genome_info then reads LastUpdateDate, RefSeq_category,
    FtpPath_RefSeq and FtpPath_GenBank off each record.
    """

    records = _live("assembly esummary",
                    lambda: _get_records([STABLE_ASSEMBLY_ID]))

    assert len(records) == 1
    record = records[0]

    # The nested DocumentSummarySet/DocumentSummary shape survived, or
    # _get_records would have raised a RuntimeError above.
    for field in ("SpeciesName",
                  "LastUpdateDate",
                  "RefSeq_category",
                  "FtpPath_RefSeq",
                  "FtpPath_GenBank"):
        assert field in record, (
            f"NCBI assembly esummary no longer has '{field}'. "
            f"topiary reads it in ncbi/entrez/proteome.py:_get_genome_info. "
            f"Fields present: {sorted(record.keys())}")

    # LastUpdateDate is slash-separated and topiary converts it to ISO before
    # parsing. If the format changed, that conversion breaks.
    raw_date = re.sub("/", "-", record["LastUpdateDate"])
    datetime.datetime.fromisoformat(raw_date)

    # RefSeq_category is compared against these two literals; anything else
    # falls through to the "other" branch, so an unknown value is not fatal --
    # but the two known values disappearing would be.
    assert isinstance(record["RefSeq_category"], str)

    # At least one ftp path has to be non-empty or topiary cannot download the
    # proteome at all.
    assert (str(record["FtpPath_RefSeq"]) != ""
            or str(record["FtpPath_GenBank"]) != ""), \
        "assembly record has neither a RefSeq nor a GenBank ftp path"
