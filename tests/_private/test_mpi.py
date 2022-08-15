
import pytest

from topiary._private.mpi import get_num_slots
from topiary._private.mpi import check_mpi_configuration

from topiary.generax import GENERAX_BINARY

import os

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_get_num_slots():

    # Should fail and catch because binary not good
    with pytest.raises(RuntimeError):
        get_num_slots("not_a_binary")

    # This will work if test environment has more than one slot
    num_slots = get_num_slots(GENERAX_BINARY)
    assert num_slots > 1



@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_check_mpi_configuration():

    # This should run fine
    check_mpi_configuration(1,GENERAX_BINARY)

    # This should run fine (assuming test environment has more than one slot)
    check_mpi_configuration(2,GENERAX_BINARY)

    # An unlikely number of slots. If this test ever passes, we've reached the
    # singularity and this code does not matter anyway
    with pytest.raises(ValueError):
        check_mpi_configuration(1000000,GENERAX_BINARY)