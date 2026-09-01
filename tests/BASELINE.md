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

From `./tests/audit_tests.py tests`, as recorded at Stage 0:

```
351 test functions: 304 checked, 29 stub, 18 assertion-free, 1 shadowed
```

**After Stage 3 this is `376 test functions: 376 checked, 0 stub, 0
assertion-free, 0 shadowed`, and `run_all_tests.sh` now runs the audit with
`--max-stub 0 --max-noassert 0`, so it is a hard gate.** The Stage 0 detail
below is kept as the record of what was wrong.

- **29 stub** — body is only `pass`. These report as passing and verify
  nothing. Includes `test_align`, `test_df_from_seed`, `test_recip_blast`,
  `test_ncbi_blast` — headline entry points of the package.
- **18 assertion-free** — the body does real work but nothing can make it
  fail. Mostly `draw/test_prettytree.py` and `reports/test_elements.py`.
- **1 shadowed** — `test_PrettyTree_render_formats` was defined twice in
  `tests/topiary/draw/test_prettytree.py`. Python kept the second; the
  shadowed first definition was the one that actually asserted
  (`assert output.exists()`). **Fixed in Stage 2**, which forced the issue: the
  tier marker landed on the shadowed copy, leaving the running copy classified
  `unit`, where the subprocess guard caught it shelling out through toyplot's
  reportlab backend. The two are now merged into one asserting `smoke` test
  covering pdf/svg/png/html. The audit now reports 0 shadowed.

Stage 3 drives stub and assertion-free toward zero. Ratchet the limits down as
that happens:

```bash
./tests/audit_tests.py tests --max-stub 0 --max-noassert 0
```

## Test tiers (Stage 2)

Every test carries exactly one tier. `integration` is derived from the
`--run-*`/`network` markers, `smoke` is explicit, `unit` is the default and is
enforced by a subprocess guard rather than trusted.

```
unit          226 tests     4.4s     hermetic
smoke          77 tests    25.2s     real data, real file I/O
integration    56 tests   185.0s     external binary or live service
```

The tiers were assigned by measurement, not by reading code: a throwaway pytest
plugin recorded per-test duration, whether `subprocess.Popen` was called, and
which data fixtures were requested, over a full all-flags run. Anything with an
opt-in marker became `integration`; of the rest, anything that spawned a
process, took over 0.25s, or pulled a real-data fixture became `smoke`.

`tests/integration-tests/` is gone — its single test was hermetic despite the
directory name and now lives at `tests/topiary/test_ete4_asr_integration.py`.

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



## Stage 3: tests that verified nothing

All 29 stubs and all 18 assertion-free tests are gone. Coverage moved with them.

```
                 Stage 0        Stage 3
lines            88.3%          92.3%    (8244/8930)
branches         84.0%          87.6%    (2991/3416)
combined         87.1%          91%
```

Worklist modules, before -> after (combined line+branch):

```
util/create_nicknames.py      12%  ->  95%
cli_scripts/seed_to_alignment      0%  -> 100%
cli_scripts/alignment_to_ancestors 0%  -> 100%
cli_scripts/create_report          0%  -> 100%
muscle/muscle.py              56%  ->  86%
io/seed.py                    48%  ->  71%
raxml/model.py                71%  ->  93%
raxml/ancestors.py            90%  ->  92%
```

How the stubs were resolved:

- **5 were redundant** -- a real mocked test already existed under a different
  name in the `*_mock.py` twin (`test_ncbi_blast`, `test__ncbi_blast_thread_function`,
  `test_merge_and_annotate`, `test__run_blast`, `test_recip_blast`). Deleted.
- **1 was an orphan** -- `tests/topiary/reports/test_quality.py` tested
  `check_duplication` in `reports/quality.py`, a module that never existed. The
  function is now `_check_duplication` in `reports/cards/duplications.py` and is
  covered by `test_duplications`. File deleted.
- **The remaining 23 got real tests**, mocked at the boundary following the
  `*_mock.py` pattern already in the repo.

Three placeholder tests (`test__thread`, `test__follow_log_generator`,
`test__follow_log_subproc_wrapper`) existed only to stop the old completeness
crawler flagging them -- they referenced the function without calling it. Since
the crawler is gone, they were replaced with tests that actually exercise the
code.

### What is deliberately NOT done

The Stage 3 plan called for folding each `*_mock.py` into its twin. **Skipped
on purpose.** The plan assumed the collisions were stub-vs-real; measuring them
showed 11 of them are real-vs-real (e.g. `test_local_blast` exists as both a
live test and a mocked one). Folding would mean renaming 11 working tests for a
purely cosmetic gain, now that the tier markers already express the
unit/integration split that the filenames were standing in for. The genuine
problem -- stubs shadowing real mock coverage -- was fixed by deleting those 5
stubs.

### Note on CLI coverage

`cli_scripts/dataframe_from_fasta.py` and `cli_scripts/fasta_into_dataframe.py`
still report 0% despite having passing tests. Those tests invoke the console
script through `subprocess.run`, so the work happens in another process and
coverage.py never sees it. The tests are real; the coverage number is an
artifact of how they run.


## Stage 4: branch coverage

```
                 Stage 0    Stage 3    Stage 4
lines            88.3%      92.3%      93.5%   (8362/8942)
branches         84.0%      87.6%      89.1%   (3046/3418)
combined         87.1%      91%        92%
```

`run_all_tests.sh` now fails if combined coverage drops below **92%**
(`coverage report --fail-under=92`). Ratchet it up as coverage improves, never
down.

Modules taken off the worklist:

```
ncbi/entrez/mrca.py       0% branch  -> 100%   (no test file existed at all)
opentree/ott.py          41% branch  ->  98%
_private/threads.py      72% branch  ->  ~95%
_private/animation.py    67% branch  ->  ~95%
reports/cards/ancestor.py 63% branch ->  100%
quality/redundancy.py    80% branch  ->  ~90%
```

### Four product bugs the missing coverage was hiding

1. **`thread_manager(shared_kwarg=...)` was broken for every documented input.**
   The type checks were written `np.issubdtype(int,to_share)` -- arguments
   backwards, so numpy tried to interpret the *value* as a dtype. Python ints,
   floats and lists all raised `TypeError`; a numpy bool array raised
   `ValueError`. It worked at all only because its single caller
   (`quality/redundancy.py`) happens to pass a numpy *integer* array, the one
   type numpy accepts as dtype-like. Fixed to check the dtype properly.

2. **Three Entrez handles were never closed** (`mrca.py`, `taxid.py`,
   `proteome.py` x2), leaking a socket per call. `sequences.py` already closed
   its handle, so the pattern existed but had not been applied consistently.
   Same class of bug as the Open Tree of Life leak.

3. **The ancestor report card rendered a literal "None"** in the "Taxonomic
   distribution of descendants" row when an ancestor had no taxonomic distance.
   The neighbouring rows already rendered "N/A", and the card header already
   guarded for None -- the table just did not.

4. **A dead assertion in `test_check_iter`.** It built a list of valid
   iterables, then looped over the *invalid* list a second time, so nothing ever
   checked that a good value is accepted. Also three `check_iter` calls sat
   inside a single `with pytest.raises` block, where only the first ever runs.

### Parametrization

`test_check_iter` was one 60-line function with 14 sequential `pytest.raises`
blocks; it is now 50 independently-reported parametrized cases. That is what
surfaced bug 4 above: sequential blocks hide both dead code and masked calls.

15 other tests still have >= 6 sequential `pytest.raises` and are worth the same
treatment (`test__prepare_for_blast` x3, `test_run_raxml`, `test_write_phy`,
`test_write_fasta`, ...).

### Still under 80% branch

These are large orchestrators; they need pipeline-level fixtures rather than
unit tests, and are the natural next target:

```
  66%   31 missed  io/seed.py
  73%   31 missed  generax/_crawl.py
  74%   18 missed  pipeline/alignment_to_ancestors.py
  74%   17 missed  pipeline/seed_to_alignment.py
```


## Stage 5: CI

The workflow is split by tier. Measured locally, running each job's exact
command:

```
unit          310 passed    9.5s   pip only -- no conda, no external binaries
smoke          89 passed     28s   conda (muscle, blast); no compile step
integration   442 passed    170s   compiled raxml-ng + generax + network
```

The `unit` job needing nothing but `pip install .[test]` is not an assumption --
it was verified by building a clean venv with the conda binaries removed from
PATH and running the whole job (flake8, the audit gate, and 310 tests) green.
Only 4 smoke tests need an external binary (3 need muscle, 1 needs toyplot's
reportlab PNG backend), which is why `smoke` takes conda but skips the compile.

### What changed and why

- **`max-parallel: 1` and the repo-wide concurrency group are gone from `unit`
  and `smoke`.** They exist now only on `integration`, which is the job that
  actually contends for NCBI and Open Tree of Life. Six matrix jobs used to run
  end to end, each compiling or restoring RAxML-NG and GeneRax, to work around
  contention that only one tier causes.
- **Only `integration` compiles the external binaries.** That was the expensive
  step in every one of the six old jobs.
- **`pip install .[dev]` fixed to `.[test]`.** `pyproject.toml` defines `test`
  and `docs`; there is no `dev` extra, so that line silently installed nothing
  and the workflow relied on a separate explicit pip install right after it.
- **`--run-network` added.** After Stage 1 the 19 network tests are opt-in, so
  without this flag CI would silently stop running them.
- **`--run-ncbi-server` is nightly-only**, to avoid NCBI throttling on every push.
- **flake8 and the audit gate moved into the fast `unit` job**, so a lint error
  or a test that verifies nothing fails in about two minutes rather than after
  the full serialized run.
- **Coverage floor enforced in the nightly job** (`--fail-under=92`), the only
  run that exercises every flag.
- **`TOPIARY_MPI_OVERSUBSCRIBE` removed.** CI set it; nothing in the repo ever
  read it. Its counterpart `TOPIARY_BLOCK_MPI_OVERSUBSCRIBE` was set and unset
  by a conftest fixture that nothing read either, guarding an empty, untracked
  `tests/topiary/_private/mpi/` directory. Both are vestiges of the MPI support
  topiary no longer has; the fixture and the directory are gone.

### Test data left in the repo, deliberately

The plan floated moving the largest `tests/data` payloads (two ~15 MB proteome
archives, a 10 MB BLAST XML) behind a download-and-cache fixture to shrink the
119 MB repo. **Not done, and it should not be.** It would put a live network
dependency back into the test suite -- the exact class of problem Stage 1 spent
its time removing, and the thing that made CI flaky in the first place. A slow
clone is a much cheaper problem than a suite that fails when a remote host is
down. If repo size becomes pressing, git-lfs keeps the inputs local to the test
run and is the better trade.


## Warnings

The suite emits **zero warnings**, and `pyproject.toml` sets
`filterwarnings = ["error", ...]` so it stays that way.

What was there before:

- **`UserWarning` from `opentree/util.py`** (5 occurrences) — "N of M species
  could not be assigned an ott id". Entirely expected: the tests deliberately
  feed unresolvable species names. Now asserted with
  `pytest.warns(UserWarning,match="could not be assigned")` instead of printed,
  which makes the warning a checked part of the contract rather than noise.

- **`RuntimeWarning: invalid value encountered in divide`** at
  `quality/alignment.py:242` and `:246` (3 occurrences). Both are 0/0 on a
  degenerate alignment: the first when every column was gaps-only and got
  dropped, the second when every column is sparse so there are no dense columns.
  NaN is the correct answer for an undefined fraction *and* is already the
  sentinel those two columns use — `score_alignment` initializes both to
  `np.nan` and fills in only the kept rows, and `polish.py`'s comparisons treat
  NaN conservatively (a NaN row is never dropped). So the division is now
  wrapped in `np.errstate(invalid="ignore")` with a comment, and
  `test_score_alignment` asserts the NaN rather than leaving it incidental.

- **A dead filter.** `"ignore:invalid escape sequence:DeprecationWarning"` was
  in `pyproject.toml` and no longer matched anything; removed.

### What `error` deliberately does not cover

`ResourceWarning` and `pytest.PytestUnraisableExceptionWarning` are ignored.
Both come from ete4 reading newick files as `data = open(data).read()`, leaving
the handle to refcounting. Under CPython the file closes immediately, so this is
style noise rather than a leak — consistent with the flat descriptor counts
measured in Stage 1. Silencing them is not hiding a real problem, and
`--fd-report` remains the tool for actual descriptor leaks. Turning them into
errors would mean rewriting every ete4 call site in topiary to read files itself,
which is churn in core tree I/O for no behavioural gain.

Note that `"error"` applies to dependency warnings too, so a numpy or pandas
deprecation on a platform not tested here will fail CI. That is the intended
trade: the fix is a one-line targeted ignore, and the warning is telling you
something real about a dependency you use.
