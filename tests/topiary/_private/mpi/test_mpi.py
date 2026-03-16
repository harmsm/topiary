
import pytest
from unittest import mock

from topiary._private.mpi import get_hosts
from topiary._private.mpi import get_num_slots
from topiary._private.mpi import check_mpi_configuration

from topiary.generax import GENERAX_BINARY

import os

@pytest.mark.run_generax
def test_get_hosts():

    hosts = get_hosts(1)
    assert len(hosts) == 1
    assert isinstance(hosts[0],str)

    hosts = get_hosts(2)
    assert len(hosts) == 2
    assert isinstance(hosts[0],str)
    assert isinstance(hosts[1],str)


@pytest.mark.run_generax
def test_get_num_slots():

    # This should work for any system with more tha one core. 
    num_slots = get_num_slots()
    assert num_slots > 1

    # Enforce a single slot with an environment variable
    with mock.patch.dict(os.environ, {'TOPIARY_MAX_SLOTS': '1'}):
        num_slots = get_num_slots()
        assert num_slots == 1


@pytest.mark.run_generax
def test_check_mpi_configuration():

    # This should run fine
    check_mpi_configuration(1)

    # This should run fine (assuming test environment has more than one slot)
    check_mpi_configuration(2)

    # An unlikely number of slots. Mock get_hosts to fail since mpirun -np 1M
    # might hang on some systems instead of failing. 
    with mock.patch("topiary._private.mpi.mpi.get_hosts",
                    side_effect=RuntimeError("fail")):
        with pytest.raises(RuntimeError):
            check_mpi_configuration(1000000)


@pytest.mark.run_generax
def test_get_num_slots_oversubscribe():

    # Reset environment to ensure no TOPIARY_MAX_SLOTS is leaking in
    if "TOPIARY_MAX_SLOTS" in os.environ:
        del os.environ["TOPIARY_MAX_SLOTS"]

    # Mock _get_mpi_oversubscribe to return True
    with mock.patch("topiary._private.mpi.mpi._get_mpi_oversubscribe", return_value=True):
        
        # Mock os.cpu_count to return 10
        with mock.patch("os.cpu_count", return_value=10):
            
            # Should return cpu count
            num_slots = get_num_slots()
            assert num_slots == 10

            # Mock TOPIARY_MAX_SLOTS to 5
            with mock.patch.dict(os.environ, {'TOPIARY_MAX_SLOTS': '5'}):
                num_slots = get_num_slots()
                assert num_slots == 5

            # Mock TOPIARY_MAX_SLOTS to 15
            with mock.patch.dict(os.environ, {'TOPIARY_MAX_SLOTS': '15'}):
                num_slots = get_num_slots()
                assert num_slots == 10
        
        # Mock os.cpu_count to return None (fallback to 1)
        with mock.patch("os.cpu_count", return_value=None):
            num_slots = get_num_slots()
            assert num_slots == 1
