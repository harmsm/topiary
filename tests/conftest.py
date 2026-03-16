import pytest
import os

def pytest_addoption(parser):
    """
    Add options to the pytest command line parser.
    """

    parser.addoption("--run-generax",
                     action="store_true",
                     default=False,
                     help="Run tests involving generax")

    parser.addoption("--run-raxml",
                     action="store_true",
                     default=False,
                     help="Run tests involving raxml")

    parser.addoption("--run-blast",
                     action="store_true",
                     default=False,
                     help="Run tests involving blast")

def pytest_collection_modifyitems(config, items):
    """
    Look for run_generax and run_raxml decorators. Modify test collection based
    on 1) pytest command line arguments and 3) operating system.
    """

    # Look for --run-generax argument. Skip test if this is not specified.
    if not config.getoption("--run-generax"):
        skipper = pytest.mark.skip(reason="Only run when --run-generax is given")
        for item in items:
            if "run_generax" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-raxml argument. Skip test if this is not specified.
    if not config.getoption("--run-raxml"):
        skipper = pytest.mark.skip(reason="Only run when --run-raxml is given")
        for item in items:
            if "run_raxml" in item.keywords:
                item.add_marker(skipper)

    # Look for --run-blast argument. Skip test if this is not specified.
    if not config.getoption("--run-blast"):
        skipper = pytest.mark.skip(reason="Only run when --run-blast is given")
        for item in items:
            if "run_blast" in item.keywords:
                item.add_marker(skipper)

    # If this is a windows box, skip any test with run_generax or run_raxml
    # decorators.
    if os.name == "nt":
        disallowed_dec = set(["run_generax","run_raxml"])
        skipper = pytest.mark.skip(reason="cannot run on windows")
        for item in items:
            if len(set(item.keywords).intersection(disallowed_dec)) > 0:
                item.add_marker(skipper)

@pytest.fixture(autouse=True)
def block_oversubscribe(request):
    """
    Ensure that TOPIARY_BLOCK_MPI_OVERSUBSCRIBE is set for heavy tasks
    (run_generax) unless they are in the mpi test directory.
    """

    # If this is a run_generax test
    if "run_generax" in request.keywords:

        # If it is NOT in the mpi test directory
        if "tests/topiary/_private/mpi" not in str(request.path):

            # Set environment variable to block oversubscribe
            os.environ["TOPIARY_BLOCK_MPI_OVERSUBSCRIBE"] = "1"
            yield
            
            # Clean up
            if "TOPIARY_BLOCK_MPI_OVERSUBSCRIBE" in os.environ:
                del os.environ["TOPIARY_BLOCK_MPI_OVERSUBSCRIBE"]
            
            return

    yield
