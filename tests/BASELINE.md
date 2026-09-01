# Test suite baseline

Recorded at the end of **Stage 0** of the test modernization work, so later
stages have something honest to measure against.

Measured from a full `run_all_tests.sh`-equivalent run (all four `--run-*`
flags enabled) on macOS, Python 3.13.15 / pytest 9.1.1 / coverage 7.15.4,
against `main` at `19b303f`.

## Why this file exists

Before Stage 0 there was no `[tool.coverage]` config, so `coverage run`
instrumented everything it imported. The reported total blended three
unrelated things:

| what was measured        | statements | covered |
| ------------------------ | ---------: | ------: |
| `src/topiary` (intended) |      8,837 |   88.3% |
| the test files themselves|      8,623 |   93.9% |
| an editable `ete4` checkout |   8,494 |   18.5% |
| **reported total**       | **25,954** | **67.3%** (badge said 63.23%) |

The badge now reports `src/topiary` only.

## Baseline numbers

Derived by re-reporting the existing full-run `.coverage` data under the new
config:

```
lines     7799 / 8833   88.3%
branches  2853 / 3396   84.0%
combined                87.1%   <- what the badge currently shows
```

**These will drop slightly on the next real full run, and that is correct.**
The old data file only contained modules that were actually imported. Now that
`source = ["src/topiary"]` is set, coverage also counts modules that are never
imported at all. Six such modules exist, all in `cli_scripts/`:

```
cli_scripts/alignment_to_ancestors.py    9 stmts   2 branches   0%
cli_scripts/bootstrap_reconcile.py       8 stmts   2 branches   0%
cli_scripts/create_report.py             9 stmts   2 branches   0%
cli_scripts/dataframe_from_fasta.py     26 stmts   4 branches   0%
cli_scripts/fasta_into_dataframe.py     13 stmts   2 branches   0%
cli_scripts/seed_to_alignment.py         9 stmts   2 branches   0%
                                        74 stmts  14 branches
```

They were invisible for exactly the reason they are untested: their tests are
`pass` stubs, so the module was never even imported. Expect the next full run
to report:

```
lines     7799 / 8907   87.6%
branches  2853 / 3410   83.7%
```

Use those as the Stage 4 starting line, not the numbers above.

Coverage is measured with `branch = true`. Branch coverage is the number
Stage 4 should move; several modules sit above 90% on lines while leaving
their error paths — the `ValueError` branches that make topiary's argument
checking trustworthy — completely unexercised.

## Weakest modules by branch coverage

Modules with at least 10 branches, worst first. This is the Stage 4 worklist.

```
  0.0%     0/14    ncbi/entrez/mrca.py
  0.0%     0/18    util/create_nicknames.py
 40.0%     4/10    raxml/tree.py
 40.9%     9/22    opentree/ott.py
 44.6%    41/92    io/seed.py
 45.0%    18/40    muscle/muscle.py
 58.7%    27/46    raxml/model.py
 70.0%    42/60    reports/reports.py
 71.7%    33/46    _private/threads.py
 73.3%    85/116   generax/_crawl.py
 73.5%    50/68    pipeline/alignment_to_ancestors.py
 74.2%    49/66    pipeline/seed_to_alignment.py
 75.0%    33/44    quality/shrink.py
 75.0%    48/64    raxml/ancestors.py
 79.5%    62/78    quality/redundancy.py
 80.0%    40/50    _private/interface.py
 81.7%    49/60    generax/_generax.py
 82.9%    68/82    _private/installed.py
```

## Test audit baseline

From `./tests/audit_tests.py tests`:

```
351 test functions: 304 checked, 29 stub, 18 assertion-free, 1 shadowed
```

- **29 stub** — body is only `pass`. These report as passing and verify
  nothing. Includes `test_align`, `test_df_from_seed`, `test_recip_blast`,
  `test_ncbi_blast` — headline entry points of the package.
- **18 assertion-free** — the body does real work but nothing can make it
  fail. Mostly `draw/test_prettytree.py` and `reports/test_elements.py`.
- **1 shadowed** — `test_PrettyTree_render_formats` is defined twice in
  `tests/topiary/draw/test_prettytree.py` (lines 291 and 406). Python keeps
  the second. The shadowed first definition is the one that actually asserts
  (`assert output.exists()`); the version that runs asserts nothing. Left
  as-is in Stage 0 because fixing it changes what executes — it belongs to
  Stage 3.

Stage 3 drives stub and assertion-free toward zero. Ratchet the limits down as
that happens:

```bash
./tests/audit_tests.py tests --max-stub 29 --max-noassert 18
```

## Hidden dependencies — status after Stage 1

- **19 of 313** default-run tests reached the live internet with no marker:
  Open Tree of Life, `ftp.ncbi.nlm.nih.gov`, and a CDN fetch inside report
  generation. **Fixed in Stage 1**: those 19 now carry `@pytest.mark.network`
  and are skipped unless `--run-network` is given, and a `block_network`
  autouse fixture blocks outbound connections for every other test so a new
  hidden dependency fails by name instead of becoming the next intermittent
  CI failure.
- **169 `os.chdir` calls** across 28 test files, with zero `try`/`finally`.
  **Fixed in Stage 1**: a `restore_cwd` autouse fixture restores the working
  directory unconditionally, and the 59 (all trailing, none mid-test) manual
  restore lines plus their `current_dir = os.getcwd()` assignments were
  removed. `tests/test_harness.py` holds regression tests that reproduce the
  original contamination and prove it is contained.

### The file-descriptor exhaustion: found and fixed

**Root cause: the `opentree` client leaks one socket per API call.**

`OTWebServiceWrapper._http_request` appends a record to `OT.ws.call_history` for
every call, and each record holds the `requests.Response` object, which pins its
socket open. Nothing ever releases them. Measured: `+1` file descriptor per Open
Tree of Life call, forever.

Across a full run the suite made a few hundred OTL calls and climbed to **226
open descriptors**. macOS's default `ulimit -n` is **256**, so the session died
partway through — under `coverage run` the extra handles tipped it over sooner,
which is why `run_all_tests.sh` failed where a bare `pytest` merely produced two
raxml failures. Individual test files never made enough calls to notice.

The fix is one line in `topiary/opentree/util.py`:

```python
OT.ws._store_api_calls = False
```

topiary never reads `call_history`, so nothing is lost. Result, at `ulimit -n 256`:

```
before:  27 tests leaked 214 descriptors, peak 226, session died
after:    5 tests leaked   9 descriptors, peak  21, 356 passed
```

**This was a product bug, not a test bug.** `seed_to_alignment` makes an OTL call
per species block, so a long pipeline run on a stock macOS box would hit the same
wall. `tests/topiary/opentree/test_ott_leak.py` guards it — one hermetic test
asserting the setting, one `network` test that fails if descriptors climb.

The remaining 9 descriptors are small, bounded, and unrelated (ftp, matplotlib,
a redundancy thread). `src/topiary/draw/ancestor_data.py:102` is the notable one:
`fig = plt.figure()` is never closed. Left for Stage 3 since the figure may be
returned to the caller.

### An earlier hypothesis that was wrong (kept as a caution)

Stage 0 recorded a guess that `mp.Manager()` in `_private/threads.py` and
`mp.Queue()` in `_private/interface.py` / `_private/animation.py` were leaking
descriptors and causing the whole-suite fd exhaustion. **Measured in Stage 1,
that is not true.** Looping each of those call sites 8-10 times gives a
perfectly flat descriptor count, with and without cleanup, including when
tracebacks are deliberately retained the way pytest retains them for failed
tests. CPython's refcounting reclaims them as soon as they go out of scope.

The cleanup was kept anyway — deterministic release is better than relying on
GC timing — but it should not be credited with fixing anything.

Adding a diagnostic instead of another guess is what actually found the real
cause above:

```bash
pytest tests --run-raxml --run-generax --run-blast --fd-report
```

`--fd-report` records descriptors opened and not closed per test and prints the
worst offenders, writing a running tally to `reports/fd-running-total.txt` as it
goes so the data survives a session that dies from descriptor exhaustion. It
never fails a test. This is what pointed at Open Tree of Life:

```
27 tests leaked 214 descriptors in total
  +56  tests/topiary/io/test_seed.py::test_read_seed
  +34  tests/topiary/opentree/test_opentree_util.py::test_ott_to_mrca
  +19  tests/topiary/pipeline/...::test_integrated_minimal_ali_to_anc
  +12  tests/topiary/quality/test_shrink.py::test_shrink_dataset
  +11  tests/topiary/opentree/test_ott.py::test_get_df_ott
  ...  (nearly every entry an Open Tree of Life caller)
```

Every large entry was an OTL caller. That is what identified the culprit — and
note it only shows up on the **full** flag set; the default suite leaked a
harmless 6.

### Scope of the network guard

The guard deliberately applies only to tests that are **not** behind an opt-in
flag. Tests marked `network`, `run_ncbi_server`, `run_blast`, `run_generax`, or
`run_raxml` are exempt.

The first cut of this guard exempted only `network`, which broke 13 tests in the
full `--run-*` run — none of which appear in the default suite, so the mistake
was invisible until someone ran `run_all_tests.sh`. The lesson is recorded here
rather than just fixed: **verify against the full flag set, not the default
suite.**

Those 13 also surfaced a genuinely surprising dependency: `run_generax` and
`run_raxml` tests contact the Open Tree of Life API, because setting up a
GeneRax run annotates the species tree from OTL. `--run-generax` therefore needs
live internet, not just the binary. That is a real constraint on CI, where the
generax job currently assumes only that the binary was compiled.

## Timings and counts

```
pytest tests                    (pre-Stage 1)   60s   318 passed,  36 skipped
pytest tests                    (post-Stage 1)  32s   299 passed,  55 skipped
pytest tests --run-network                      60s   318 passed,  36 skipped
pytest tests <all five flags>                  178s   356 passed,   0 skipped
```

The default suite got roughly twice as fast purely by not waiting on remote
servers, and nothing was lost — `--run-network` reproduces the old result
exactly. The all-flags run is the number `run_all_tests.sh` produces and the one
to check before claiming a change is safe. Check it at a realistic descriptor
limit too -- `bash -c 'ulimit -n 256; pytest ...'` -- since a developer box
defaults to 256 and a leak that is invisible at a high limit is fatal there.

