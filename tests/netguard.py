"""
Support types for the network guard in tests/conftest.py.

These live here rather than in conftest.py so they can be pickled. There are two
conftest.py files in this repo (tests/ and tests/topiary/), both of which pytest
imports under the module name "conftest". An exception class defined in one of
them therefore has an ambiguous __module__, and multiprocessing cannot send it
back to the parent process:

    PicklingError: Can't pickle <class 'conftest.NetworkAccessBlocked'>:
    attribute lookup NetworkAccessBlocked on conftest failed

That matters because the guard *can* fire inside a worker process: on Linux the
default multiprocessing start method is fork, so a child inherits the patched
socket.connect. Without a picklable exception, a test that genuinely touches the
network from a pool worker fails with an unreadable PicklingError instead of a
message naming the test and the address.
"""

import socket


class NetworkAccessBlocked(Exception):
    """
    A test that is not marked `network` tried to open an outbound connection.
    """

    pass


# Addresses a test may always reach: loopback, and anything on the local
# machine. multiprocessing uses AF_UNIX sockets for its Manager and Queue
# plumbing, so those must stay open or thread_manager tests cannot run.
LOCAL_HOSTS = frozenset(["127.0.0.1", "::1", "localhost", "0.0.0.0", ""])


def is_local_address(sock, address):
    """
    True if `address` is loopback or a local unix socket.
    """

    if getattr(sock, "family", None) == socket.AF_UNIX:
        return True

    # AF_UNIX addresses are plain paths
    if isinstance(address, (str, bytes)):
        return True

    if isinstance(address, (tuple, list)) and len(address) > 0:
        return address[0] in LOCAL_HOSTS

    return False
