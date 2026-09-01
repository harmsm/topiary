
import pytest

from topiary._private import threads
from topiary._private.threads import _thread
import time, random

import numpy as np

def test_MockLock():

    # MockLock is substituted for a real multiprocessing lock when we drop to a
    # single thread, so it has to satisfy the same interface without doing
    # anything. The acquire() return value matters: callers check it.
    lock = threads.MockLock()

    assert lock.acquire() is True
    assert lock.acquire(blocking=False) is True
    assert lock.acquire(blocking=True,timeout=0.1) is True
    assert lock.release() is None

    # Has to work as a plain acquire/release pair the way calling code uses it
    lock.acquire()
    lock.release()

def test_get_num_threads():

    num_threads = threads.get_num_threads(-1,manual_num_cores=10)
    assert num_threads == 10

    num_threads = threads.get_num_threads(1,manual_num_cores=10)
    assert num_threads == 1

    num_threads = threads.get_num_threads(5,manual_num_cores=10)
    assert num_threads == 5

    num_threads = threads.get_num_threads(5,manual_num_cores=2)
    assert num_threads == 2

    bad_num_threads = [0,-2,"test",None,int,1.1]
    for b in bad_num_threads:
        print(f"testing {b} in get_num_threads")
        with pytest.raises(ValueError):
            threads.get_num_threads(b,manual_num_cores=10)

def _function_to_thread(value,some_option=False,lock=None):
    """
    Function for testing thread_manager.
    """

    if lock is not None:

        # If lock is acquired, the output value will change. (Not directly
        # testing locking -- that will depend too much on specific function).
        lock.acquire()
        try:
            value = 2*value
        finally:
            lock.release()

    if some_option:
        value = 3*value

    # Sleep for a random time to knock multithreaded-order out of sync
    time.sleep(random.choice([0,0.01]))

    return value

@pytest.mark.smoke
def test_thread_manager():

    kwargs_list = [{"value":2,"some_option":False},
                   {"value":2,"some_option":True},
                   {"value":3,"some_option":False}]

    # Run on single thread. (This will not use pool)
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=False)
    assert np.array_equal(output,[2,6,3])

    # Run on single thread. (This will not use pool; make sure progress bar at
    # least accepted).
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=True)
    assert np.array_equal(output,[2,6,3])

    # Run twice on two threads to make sure order is preserved
    for i in range(2):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=True,
                                        pass_lock=False)
        assert np.array_equal(output,[2,6,3])

    # Run twice on two threads make sure order is preserved
    for i in range(2):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=False,
                                        pass_lock=True)
        assert np.array_equal(output,[4,12,6])


def test__thread():
    """
    _thread runs one function call and puts (calc_number, result) on the queue.
    thread_manager relies on that calc_number to put results back in order.
    """

    import queue as queue_module

    q = queue_module.Queue()

    def double(x):
        return x*2

    _thread((3,double,{"x":21},q))

    assert q.qsize() == 1
    calc_number, result = q.get()
    assert calc_number == 3
    assert result == 42

    # A second call adds to the queue rather than replacing
    _thread((0,double,{"x":1},q))
    assert q.qsize() == 1
    assert q.get() == (0,2)


def _sum_shared(x,shared=None,lock=None):
    """
    Read a shared value/array and combine it with x. Used to check that
    shared_kwarg actually reaches the worker.
    """

    if hasattr(shared,"__iter__") or hasattr(shared,"__getitem__"):
        try:
            return x + sum(shared)
        except TypeError:
            return x + shared.value

    return x + shared.value


@pytest.mark.smoke
def test_thread_manager_shared_int():
    """
    shared_kwarg with a scalar int gets wrapped in a manager Value and handed to
    every worker.
    """

    kwargs_list = [{"x":1,"shared":10},{"x":2,"shared":10},{"x":3,"shared":10}]

    output = threads.thread_manager(kwargs_list,
                                    _sum_shared,
                                    num_threads=2,
                                    progress_bar=False,
                                    shared_kwarg="shared")

    assert np.array_equal(output,[11,12,13])


@pytest.mark.smoke
def test_thread_manager_shared_float():

    kwargs_list = [{"x":1.0,"shared":0.5},{"x":2.0,"shared":0.5}]

    output = threads.thread_manager(kwargs_list,
                                    _sum_shared,
                                    num_threads=2,
                                    progress_bar=False,
                                    shared_kwarg="shared")

    assert np.allclose(output,[1.5,2.5])


@pytest.mark.smoke
def test_thread_manager_shared_iterable():
    """
    An iterable shared_kwarg becomes a manager Array rather than a Value.
    """

    kwargs_list = [{"x":1,"shared":[1,2,3]},{"x":2,"shared":[1,2,3]}]

    output = threads.thread_manager(kwargs_list,
                                    _sum_shared,
                                    num_threads=2,
                                    progress_bar=False,
                                    shared_kwarg="shared")

    assert np.array_equal(output,[7,8])


@pytest.mark.smoke
def test_thread_manager_shared_kwarg_validation():

    # shared_kwarg names something that is not in kwargs_list
    with pytest.raises(ValueError):
        threads.thread_manager([{"x":1}],
                               _sum_shared,
                               num_threads=2,
                               progress_bar=False,
                               shared_kwarg="not_a_key")

    # An iterable of the wrong type
    with pytest.raises(ValueError):
        threads.thread_manager([{"x":1,"shared":["a","b"]}],
                               _sum_shared,
                               num_threads=2,
                               progress_bar=False,
                               shared_kwarg="shared")

    # A scalar of the wrong type
    with pytest.raises(ValueError):
        threads.thread_manager([{"x":1,"shared":"a_string"}],
                               _sum_shared,
                               num_threads=2,
                               progress_bar=False,
                               shared_kwarg="shared")


def test_get_num_threads_cpu_count_fallbacks(mocker):
    """
    os.sched_getaffinity is Linux-only and mp.cpu_count can raise, so there is a
    fallback chain down to os.cpu_count and finally to a single thread.
    """

    # sched_getaffinity missing -> falls through to mp.cpu_count
    mocker.patch.object(threads.os,"sched_getaffinity",
                        side_effect=AttributeError,create=True)
    mocker.patch.object(threads.mp,"cpu_count",return_value=8)
    assert threads.get_num_threads(-1) == 8

    # mp.cpu_count not implemented -> falls through to os.cpu_count
    mocker.patch.object(threads.mp,"cpu_count",side_effect=NotImplementedError)
    mocker.patch.object(threads.os,"cpu_count",return_value=4)
    assert threads.get_num_threads(-1) == 4

    # Nothing could work it out -> single thread
    mocker.patch.object(threads.os,"cpu_count",return_value=None)
    assert threads.get_num_threads(-1) == 1
