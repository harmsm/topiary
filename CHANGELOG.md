# Changelog

This changelog starts at v1.1.0. For anything earlier, see the git history.

## v1.1.0

39 commits since v1.0.0. The headline changes are a rewritten bootstrap
reconciliation that scales across cluster nodes without MPI, support for
microbial datasets, a new UniProt seed importer, and a substantially rebuilt
test suite that surfaced several real bugs along the way.

### Highlights

- **Bootstrap reconciliation no longer uses MPI.** `topiary-bootstrap-reconcile`
  is now a re-entrant filesystem crawler: launch as many copies as you like
  against the same run directory (a SLURM job array with one task per node, for
  example) and they coordinate through marker files. Re-running resumes rather
  than restarting.
- **Microbial datasets are supported.** `read_seed` decides automatically
  whether a dataset should be treated as species-aware, so bacterial seeds no
  longer fail on unresolvable Open Tree of Life lookups.
- **New `topiary-uniprot-to-seed`** converts a UniProt FASTA download directly
  into a topiary seed dataframe.
- **Better failure messages for broken external binaries.** A RAxML-NG or
  GeneRax binary compiled for a different CPU now reports the actual signal
  (`killed by signal SIGILL`) and says the architecture is likely wrong, instead
  of implying a `$PATH` problem.
- **A long pipeline run no longer exhausts file descriptors.** See *Bug fixes*.
- **`import topiary` went from 1.28s to 0.40s** — `topiary.draw` (matplotlib,
  toytree, ete4, and the scipy that ete4 drags in) is now loaded on first
  attribute access. `topiary.draw.tree(...)` works exactly as before.

### Breaking changes

**`bootstrap_reconcile` / `topiary-bootstrap-reconcile`** — the thread arguments
are gone, replaced by a launcher prefix. GeneRax parallelism is MPI-only, and
you now own the launcher and the allocation:

```
removed:  --num_threads  --threads_per_replicate  --restart  --overwrite
added:    --generax_launch          e.g. "mpirun -np 8"; "" runs single-core
          --replicate_timeout_factor
          --replicate_max_hours
          --replicate_min_seconds
```

`--restart` and `--overwrite` are gone because the crawler is always resumable:
re-run it to continue, and delete the `*_reconciled-tree-bootstraps` directory
to start over.

**`generax.reconcile()`** — `num_threads`, `threads_per_rep` and
`converge_cutoff` removed; `generax_launch` added.

**MPI support removed from topiary itself.** `topiary._private.mpi` and
`installed.check_mpirun()` are gone, and `mpi4py` is no longer a dependency.
`openmpi` remains in `environment.yml` only because the bundled GeneRax binary
links against it.

**Python is pinned to `<3.14`**, to avoid conda resolution failures on Linux.

### New features

- `topiary.io.uniprot_to_seed()` and the `topiary-uniprot-to-seed` CLI script.
- `read_seed(species_aware=...)` — `True` requires every species to resolve on
  the Open Tree of Life synthetic tree, `False` skips resolution entirely, and
  the default `None` picks automatically (microbial datasets get `False`).
  Threaded through `seed_to_alignment`.
- `seed_to_alignment(blast_taxid=...)` — limit the NCBI BLAST search to a taxid
  explicitly, rather than always inferring it from the key species.
- `topiary.ncbi.entrez.get_mrca_taxid()` — most recent common ancestor taxid for
  a list of species.
- `install.sh` gained `--no-generax`, `--no-raxml` and `--no-cluster`, and was
  reworked so it installs cleanly in Google Colab and on clusters.
- `alignment_to_ancestors(generax_launch=...)`, matching the new
  reconciliation interface.

### Bug fixes

- **Open Tree of Life calls leaked a socket each.** The `opentree` client
  records every API response in `OT.ws.call_history`, and each record pins the
  underlying `requests.Response` — and therefore its socket — open forever.
  Nothing released them. A `seed_to_alignment` run makes an OTL call per species
  block, so on a machine with the stock macOS limit of 256 file descriptors a
  long run would die with `OSError: [Errno 24] Too many open files`. topiary now
  turns that recording off. **This is the most likely bug to have affected a
  real analysis.**
- **Three Entrez handles were never closed** (`ncbi/entrez/mrca.py`,
  `taxid.py`, and two in `proteome.py`), leaking a socket per call in the same
  way.
- **`align()` left scratch files behind on failure.** The temporary
  `topiary-tmp_*_align-in.fasta` was deleted only after muscle returned, so a
  muscle crash stranded it in the caller's working directory. Cleanup now
  happens in a `finally`.
- **`thread_manager(shared_kwarg=...)` was broken for every documented input.**
  The type checks were written `np.issubdtype(int, value)` — arguments reversed,
  so numpy tried to interpret the *value* as a dtype. Python ints, floats and
  lists all raised `TypeError`, and a numpy bool array raised `ValueError`. It
  worked at all only because its single caller happens to pass a numpy integer
  array, the one type numpy accepts as dtype-like.
- **The ancestor report card printed a literal `None`** in the "Taxonomic
  distribution of descendants" row for an ancestor with no taxonomic
  distribution, where the neighbouring rows already rendered `N/A`.
- **Report generation crashed on dataframes with no `ott` column**
  (`reports/cards/input.py` checked the column's contents before checking it
  existed).
- **`MockLock.acquire()` did not accept `blocking`/`timeout`**, so it was not a
  drop-in for a real `multiprocessing` lock.
- **`multiprocessing` resources are now released deterministically** — the
  `Manager` in `threads.py` and the `Queue`/`Process` pairs in `interface.py`
  and `animation.py` are shut down explicitly rather than left to the garbage
  collector.

### Testing and CI

Not user-facing, but it is most of the diff and it is why several of the fixes
above were found.

- **417 test functions (451 collected cases), up from 318, and every one of
  them verifies something.** At v1.0.0 the suite contained 29 `pass`-bodied
  stubs, 15 tests with no assertion at all, and one test silently shadowed by a
  duplicate function name — 45 of 318 reporting green while checking nothing.
  All are now real tests, and `tests/audit_tests.py` fails the build if any
  reappear.
- **Coverage is measured correctly for the first time.** The badge previously
  blended `src/topiary` with an installed `ete4` checkout and the test files
  themselves. It now measures topiary only: **93.5% of lines, 89.1% of
  branches**, with a 92% floor enforced in CI.
- **Tests are tiered** — `pytest -m unit` runs 311 hermetic tests in about ten
  seconds, `-m smoke` adds real file I/O, and `integration` covers the external
  binaries and live services. The `unit` tier is enforced rather than trusted: a
  unit test that shells out or opens a socket fails.
- **19 tests were silently depending on live internet** with no marker. They are
  now explicit and opt-in via `--run-network`, and a guard blocks outbound
  connections everywhere else, so a new hidden dependency fails by name instead
  of becoming an intermittent CI failure.
- **CI runs in parallel again.** The whole workflow used to be serialized to
  work around contention that only the integration tier causes; `unit` and
  `smoke` now run the full matrix concurrently, and only `integration` compiles
  RAxML-NG and GeneRax. It also validates a cached binary before trusting it and
  rebuilds if it will not run on that runner — the likely cause of the
  intermittent SIGILL failures.
- 108 RAxML/GeneRax log fixtures that a blanket `*.log` rule in `.gitignore` had
  been silently excluding are now committed, so tests depending on them no
  longer pass locally and fail on a fresh checkout.
