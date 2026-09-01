#!/usr/bin/env python
"""
Audit the test suite for tests that report as passing without verifying
anything.

This replaces the old completeness_crawler.py, which tried to infer whether a
source function was tested by matching test names against source names
(topiary.a.b.fcn -> tests.topiary.a.test_b.test_fcn). That heuristic was both
noisy and redundant: branch coverage answers "is this code exercised?" directly
and correctly, and the name-matching broke as soon as tests stopped being one
file per module.

What coverage cannot see is a test that runs but checks nothing -- a `pass`
body executes just fine, and a test with no assert happily "passes" whatever
the code does. That is the one real signal worth keeping, so it is what this
script reports:

    STUB      the body is only `pass` / `...` / a docstring
    NOASSERT  the body does something, but nothing can make it fail

Usage
-----
    ./tests/audit_tests.py tests/
    ./tests/audit_tests.py tests/ --max-stub 0 --max-noassert 5

With --max-stub / --max-noassert the script exits non-zero when the count
exceeds the limit, so it can be used as a CI gate that ratchets downward.
"""

import argparse
import ast
import os
import sys

# Calls that can make a test fail, and so count as a real check.
FAILING_CALLS = {
    "raises", "warns", "deprecated_call", "fail", "xfail", "exit",
    "approx", "raises_group",
    # pytester's result object: these raise on mismatch, so a test using them
    # can genuinely fail even with no bare `assert` in sight.
    "fnmatch_lines", "re_match_lines", "no_fnmatch_line",
}

# Any attribute call whose name starts with one of these counts as an
# assertion: mock's assert_called_with etc., and the numpy/pandas testing
# helpers (np.testing.assert_allclose, pd.testing.assert_frame_equal, ...).
ASSERT_PREFIXES = ("assert",)


def _call_name(node):
    """
    Return the dotted-ish name of a Call's func, e.g. "pytest.raises" -> both
    "raises" (attr) and "pytest.raises" (full) as a 2-tuple.
    """

    func = node.func
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id

    return None


def _is_trivial(stmt):
    """
    True if a statement does nothing: `pass`, a bare `...`, or a docstring.
    """

    if isinstance(stmt, ast.Pass):
        return True

    if isinstance(stmt, ast.Expr):
        value = stmt.value
        if isinstance(value, ast.Constant):
            # Covers both a docstring and a bare `...`
            return True

    return False


def checking_helpers(tree):
    """
    Names of module-level helper functions whose bodies can fail.

    A test that calls `_load_check_html(out)` is checked, even though the
    assertion lives one frame down. Without this, factoring a check out into a
    helper -- which is good practice -- would look like a missing assertion.
    Only one level deep, which covers the pattern in this codebase.
    """

    helpers = set()

    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if node.name.startswith("test_"):
            continue
        if _body_can_fail(node):
            helpers.add(node.name)

    return helpers


def _body_can_fail(node):
    """
    True if anything in this function's body could raise an assertion-style
    failure.
    """

    for sub in ast.walk(node):

        if isinstance(sub, ast.Assert):
            return True

        if isinstance(sub, ast.Raise):
            return True

        if isinstance(sub, ast.Call):
            name = _call_name(sub)
            if name is None:
                continue
            if name in FAILING_CALLS or name.startswith(ASSERT_PREFIXES):
                return True

    return False


def classify(node, helpers=frozenset()):
    """
    Classify a test function node as "stub", "noassert", or "ok".

    Parameters
    ----------
    node : ast.FunctionDef or ast.AsyncFunctionDef
        the test function to classify
    helpers : set of str
        names of module-level helpers that can themselves fail; calling one
        counts as checking

    Returns
    -------
    verdict : str
        one of "stub", "noassert", "ok"
    """

    body = [s for s in node.body if not _is_trivial(s)]
    if len(body) == 0:
        return "stub"

    for sub in ast.walk(node):

        # A plain `assert something`
        if isinstance(sub, ast.Assert):
            return "ok"

        # pytest.raises(...), pytest.fail(...), mock.assert_called_once(),
        # np.testing.assert_allclose(), ...
        if isinstance(sub, ast.Call):
            name = _call_name(sub)
            if name is None:
                continue
            if name in FAILING_CALLS:
                return "ok"
            if name.startswith(ASSERT_PREFIXES):
                return "ok"
            if name in helpers:
                return "ok"

        # `with pytest.raises(...)` used without a call we caught above, and
        # bare `raise` inside a test body (deliberate failure signalling).
        if isinstance(sub, ast.Raise):
            return "ok"

    return "noassert"


def collect(tree):
    """
    Yield the test functions pytest would actually collect from a parsed
    module: functions named test_* at module level, and methods named test_*
    on classes named Test*.

    Functions nested inside other functions are deliberately skipped. pytest
    does not collect them, and they are usually little helper callables that a
    real test defines and then asserts on -- counting them produces false
    "assertion-free test" reports.
    """

    for node in tree.body:

        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("test_"):
                yield node

        elif isinstance(node, ast.ClassDef) and node.name.startswith("Test"):
            for sub in node.body:
                if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    if sub.name.startswith("test_"):
                        yield sub


def audit(test_dir):
    """
    Walk test_dir and classify every test function found.

    Parameters
    ----------
    test_dir : str
        directory to crawl

    Returns
    -------
    results : dict
        keys "stub", "noassert", "ok"; values are lists of
        "path:lineno name" strings
    """

    results = {"stub": [], "noassert": [], "ok": [], "duplicate": []}

    for dirpath, dirs, files in os.walk(test_dir):

        # Don't descend into caches or virtualenv-ish directories
        dirs[:] = [d for d in dirs if d != "__pycache__" and not d.startswith(".")]

        for f in sorted(files):

            if not (f.startswith("test_") and f.endswith(".py")):
                continue

            path = os.path.join(dirpath, f)

            try:
                tree = ast.parse(open(path).read(), filename=path)
            except SyntaxError as e:
                print(f"SYNTAX-ERROR {path}: {e}", file=sys.stderr)
                continue

            helpers = checking_helpers(tree)

            seen = {}
            for node in collect(tree):

                # A second def with the same name silently replaces the first,
                # so the earlier test never runs. pytest reports no error --
                # it just collects one fewer test than the file appears to
                # contain.
                if node.name in seen:
                    results["duplicate"].append(
                        f"{path}:{node.lineno} {node.name} "
                        f"(shadows definition at line {seen[node.name]})")
                seen[node.name] = node.lineno

                verdict = classify(node, helpers)
                results[verdict].append(f"{path}:{node.lineno} {node.name}")

    return results


def main(argv=None):

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("test_dir", nargs="?", default="tests",
                        help="directory to crawl (default: tests)")
    parser.add_argument("--max-stub", type=int, default=None,
                        help="exit non-zero if more than this many stub tests are found")
    parser.add_argument("--max-noassert", type=int, default=None,
                        help="exit non-zero if more than this many assertion-free tests are found")
    parser.add_argument("--quiet", action="store_true",
                        help="print the summary only, not the individual tests")

    args = parser.parse_args(argv)

    results = audit(args.test_dir)

    n_stub = len(results["stub"])
    n_noassert = len(results["noassert"])
    n_ok = len(results["ok"])
    n_dupe = len(results["duplicate"])
    n_total = n_stub + n_noassert + n_ok

    if not args.quiet:
        for entry in results["stub"]:
            print(f"STUB       {entry}")
        for entry in results["noassert"]:
            print(f"NOASSERT   {entry}")
        for entry in results["duplicate"]:
            print(f"DUPLICATE  {entry}")
        if n_stub or n_noassert or n_dupe:
            print()

    print(f"{n_total} test functions: "
          f"{n_ok} checked, {n_stub} stub, {n_noassert} assertion-free, "
          f"{n_dupe} shadowed by a duplicate name")

    status = 0
    if args.max_stub is not None and n_stub > args.max_stub:
        print(f"FAIL: {n_stub} stub tests exceeds limit of {args.max_stub}")
        status = 1
    if args.max_noassert is not None and n_noassert > args.max_noassert:
        print(f"FAIL: {n_noassert} assertion-free tests exceeds limit of {args.max_noassert}")
        status = 1
    if n_dupe > 0:
        # Always an error: a shadowed test is silently not running.
        print(f"FAIL: {n_dupe} test(s) shadowed by a duplicate name")
        status = 1

    return status


if __name__ == "__main__":
    sys.exit(main())
