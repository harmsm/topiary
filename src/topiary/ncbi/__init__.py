"""
Interface to NCBI blast and entrez.
"""

from Bio import Entrez
Entrez.email = "topiary.phylogenetics@gmail.com"

import os
import threading
import time

# Figure out how often we can hit the NCBI servers with requests. Limit is 
# 3 times per second without an api key; 10 times per second with an api key. 
# NCBI_REQUEST_FREQ coordinates between threads so we do not exceed this limit 
try:
    Entrez.api_key = os.environ['NCBI_API_KEY']
    NCBI_REQUEST_FREQ = 1/9.5
except KeyError:
    NCBI_REQUEST_FREQ = 1/2.5


# Time of the most recent NCBI request, and a lock guarding it.
_last_request_time = 0.0
_rate_limit_lock = threading.Lock()


def rate_limit():
    """
    Block until enough time has passed to make another NCBI request.

    NCBI returns HTTP 429 ("Too Many Requests") if it is queried faster than
    3 times per second (10 with an API key). Call this immediately before every
    Entrez request.

    Unlike a bare ``time.sleep(NCBI_REQUEST_FREQ)``, this only waits for the
    remaining part of the interval, so a caller that was already slow pays
    nothing.

    Notes
    -----
    This coordinates threads within one process. Code that fans Entrez requests
    out over a multiprocessing pool (``entrez.get_sequences``,
    ``blast.ncbi_blast``) additionally passes a ``multiprocessing`` lock down to
    its workers and does its own sleep; those paths are already limited and do
    not call this function.
    """

    global _last_request_time

    with _rate_limit_lock:

        wait = NCBI_REQUEST_FREQ - (time.monotonic() - _last_request_time)
        if wait > 0:
            time.sleep(wait)

        _last_request_time = time.monotonic()

from ._parse_ncbi_line import parse_ncbi_line
from .blast import local_blast, ncbi_blast, recip_blast, make_blast_db
from .blast import records_to_df, read_blast_xml
from .blast import merge_blast_df, merge_and_annotate
from .entrez import get_sequences, get_taxid, get_proteome, get_mrca_taxid
