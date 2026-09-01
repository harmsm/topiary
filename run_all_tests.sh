#!/bin/bash

echo "Running flake8"
flake_test=`flake8 src/topiary tests --count --select=E9,F63,F7,F82 --show-source --statistics`
if [[ "${flake_test}" != 0 ]]; then
    echo "flake failed"
    exit
fi

rm -rf reports
mkdir reports
mkdir reports/junit
mkdir reports/coverage
mkdir reports/badges

echo "Running flake8, aggressive"
flake8 src/topiary tests --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics > reports/flake.txt

echo "Checking for test data that exists locally but is not in git"
# A test that reads a file which is on disk here but not committed passes
# locally and fails on a fresh checkout -- i.e. in CI. This has happened
# (a blanket *.log rule in .gitignore hid 108 raxml log fixtures), so check.
untracked_data=$(comm -23 \
    <(find tests/data -type f | sort) \
    <(git ls-files tests/data | sort) \
    | grep -vE '(\.DS_Store|\.ipynb_checkpoints/)' || true)
if [[ -n "${untracked_data}" ]]; then
    echo "FAIL: test data on disk but not tracked by git."
    echo "Tests reading these pass here and fail on a fresh checkout:"
    echo "${untracked_data}"
    exit 1
fi

echo "Auditing tests that verify nothing"
# Reports stub (`pass`-bodied) tests, tests with no assertion, and tests
# silently shadowed by a duplicate function name. As of Stage 3 all three are
# at zero, so the limits are set to zero and this is a hard gate: a new test
# that verifies nothing fails the build.
./tests/audit_tests.py tests --max-stub 0 --max-noassert 0 > reports/test-audit.txt
tail -1 reports/test-audit.txt
grep DUPLICATE reports/test-audit.txt

echo "Running coverage.py"
# Note: `source`, `branch`, and the html/xml output paths all come from
# [tool.coverage.*] in pyproject.toml, so they are not repeated here.
coverage erase
coverage run -m pytest tests --run-network --run-raxml --run-generax --run-blast --run-ncbi-server --junit-xml=reports/junit/junit.xml

echo "Generating reports"
coverage report > reports/coverage/coverage.txt

# Coverage floor. Set to the measured value at the end of Stage 4; ratchet it
# up as coverage improves, never down. This is a gate, not a badge: a change
# that drops coverage below the floor fails here.
coverage report --fail-under=92 > /dev/null
if [[ $? -ne 0 ]]; then
    echo "FAIL: coverage dropped below the floor of 92%"
    coverage report | tail -2
    exit 1
fi
coverage html
coverage xml

genbadge tests -i reports/junit/junit.xml -o docs/badges/tests-badge.svg
genbadge coverage -i reports/coverage/coverage.xml -o docs/badges/coverage-badge.svg

wget https://github.com/harmslab/topiary/actions/workflows/python-app.yml/badge.svg -O docs/badges/ghwf.svg
wget https://readthedocs.org/projects/topiary-asr/badge/?version=latest -O docs/badges/rtd.svg
