"""
Functions for interfacing with and validating mpi configuration.
"""

import subprocess
import sys
import os

from topiary._private.environment import load_env_variable
from topiary._private.check import check_bool
from topiary._private.check import check_int

def _get_mpi_oversubscribe():
    """
    Check if we should use --oversubscribe with mpirun.
    """
    v = load_env_variable("TOPIARY_MPI_OVERSUBSCRIBE",
                          check_function=check_int,
                          check_function_kwargs={"minimum_allowed":0,
                                                 "maximum_allowed":1})
    if v is None:
        return False
    
    # If we are blocking oversubscribe (set by test suite for heavy tasks)
    # return False. 
    block = load_env_variable("TOPIARY_BLOCK_MPI_OVERSUBSCRIBE",
                              check_function=check_int,
                              check_function_kwargs={"minimum_allowed":0,
                                                     "maximum_allowed":1})
    if block:
        return False

    return bool(v)

def get_mpi_env():
    """
    Get a copy of the current os.environ stripped of variables that can 
    conflict with mpirun. This is particularly important for SLURM HPC 
    environments where multiple mpirun instances launched via python 
    multiprocessing can collide.

    Returns
    -------
    env : dict
        environment dictionary suitable for passing to subprocess.run
    """
    return os.environ.copy()

def get_hosts(num_slots):
    """
    Get the hosts allocated by the job manager by running a script with mpirun
    that returns the host names.

    Returns
    -------
    host : list
        list of host names (for example: ["n001","n001","n002"])
    """

    location = os.path.dirname(os.path.realpath(__file__))
    script = os.path.join(location,"_get_hosts.py")
    python = sys.executable

    cmd = ["mpirun"]
    if _get_mpi_oversubscribe():
        cmd.append("--oversubscribe")
    cmd.extend(["-np",f"{num_slots}",python,script])

    ret = subprocess.run(cmd, capture_output=True, env=get_mpi_env())
    if ret.returncode == 0:
        stdout = ret.stdout.decode()
        raw_hosts = [s.strip() for s in stdout.split("\n") if s.strip() != ""]
        raw_hosts.sort()

        # Check if the HPC's hostname is returned. If we are running on a 
        # single-node allocation, we map the local hostname to "localhost" 
        # to prevent OpenMPI from attempting to SSH to itself, which is 
        # often explicitly blocked on compute nodes and crashes.
        import socket
        local_hostname = socket.gethostname().split(".")[0]
        hosts = []
        for h in raw_hosts:
            if h.split(".")[0] == local_hostname or h == "localhost":
                hosts.append("localhost")
            else:
                hosts.append(h)
    else:
        err = "Could not determine hosts. _get_hosts.py script returned:\n\n"
        err += f"stdout:\n\n{ret.stdout.decode()}\n\n"
        err += f"stderr:\n\n{ret.stderr.decode()}\n\n"
        raise RuntimeError(err)

    return hosts

def get_num_slots():
    """
    Get the number of mpi slots available by running get_hosts until it throws
    an error. 

    Returns
    -------
    num_slots : int
        number of slots available
    """

    # Get environment variable if defined
    max_num_slots = load_env_variable("TOPIARY_MAX_SLOTS",
                                      check_function=check_int,
                                      check_function_kwargs={"minimum_allowed":1})
    
    # If oversubscribe is enabled, return the number of cores on the machine. 
    # This prevents an infinite loop because oversubscribe will almost always 
    # work regardless of the number of slots requested. 
    if _get_mpi_oversubscribe():
        num_slots = os.cpu_count()
        if num_slots is None:
            num_slots = 1
            
        if max_num_slots is not None:
            if num_slots > max_num_slots:
                num_slots = max_num_slots
        
        return num_slots

    # Increase number of slots until mpirun fails
    num_slots = 1
    while True:

        try:
            _ = get_hosts(num_slots)
            num_slots += 1
        except RuntimeError as error:
            num_slots = num_slots - 1

            # If we have no slots, something went terribly wrong.
            if num_slots < 1:
                err = "\nCould not determine the number of MPI slots. The test script\n"
                err += "raised the following error."
                raise RuntimeError(err) from error

            break

        # Hard cap based on environment variable
        if max_num_slots is not None:
            if num_slots >= max_num_slots:
                num_slots = max_num_slots
                break


    return num_slots


def check_mpi_configuration(num_threads):
    """
    Make sure mpi configuration allows the requested number of threads.

    Parameters
    ----------
    num_threads : int
        number of threads (e.g. slots) to test. if -1, try to infer the number
        of slots using get_num_slots
    """

    # if threads were not passed in directly, infer from the environment
    if num_threads == -1:
        num_threads = get_num_slots()

    try:
        get_hosts(num_threads)
    except RuntimeError as error:
        err = "\n\nmpirun is not working. This could because you\n"
        err += "set num_threads to be more than the number of nodes you have\n"
        err += "allocated. If you did not set num_threads specifically, try\n"
        err += "setting it rather than having topiary try to figure out the\n"
        err += "number of processors. Another issue could be subtle problems\n"
        err += "with how processors are being requested via your job management\n"
        err += "software (i.e. SLURM, TORQUE, etc.). Maybe play with flags like\n"
        err += "--ntasks-per-node or talk to your cluster administrator.\n"

        raise RuntimeError(err) from error
