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

echo "Auditing tests that verify nothing"
# Reports stub (`pass`-bodied) tests, tests with no assertion, and tests
# silently shadowed by a duplicate function name. Exits non-zero only on
# duplicates for now; --max-stub/--max-noassert can be ratcheted down as
# these get fixed.
./tests/audit_tests.py tests > reports/test-audit.txt
tail -1 reports/test-audit.txt
grep DUPLICATE reports/test-audit.txt

echo "Running coverage.py"
# Note: `source`, `branch`, and the html/xml output paths all come from
# [tool.coverage.*] in pyproject.toml, so they are not repeated here.
coverage erase
coverage run -m pytest tests --run-network --run-raxml --run-generax --run-blast --run-ncbi-server --junit-xml=reports/junit/junit.xml

echo "Generating reports"
coverage report > reports/coverage/coverage.txt
coverage html
coverage xml

genbadge tests -i reports/junit/junit.xml -o docs/badges/tests-badge.svg
genbadge coverage -i reports/coverage/coverage.xml -o docs/badges/coverage-badge.svg

wget https://github.com/harmslab/topiary/actions/workflows/python-app.yml/badge.svg -O docs/badges/ghwf.svg
wget https://readthedocs.org/projects/topiary-asr/badge/?version=latest -O docs/badges/rtd.svg
