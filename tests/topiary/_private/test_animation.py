
import pytest

from topiary._private import animation

import time

@pytest.mark.smoke
def test_WaitingAnimation():

    # Not really testing animation output, but am testing threading. This test
    # really only makes sense if we're using pytest-timeout with the
    # --timeout flag. If threading doesn't start properly, time.sleep will not
    # occur; of w.stop() does not exit, this function will hang forever. 

    w = animation.WaitingAnimation(delay=0.1)

    w.start()
    assert w._proc is not None
    assert w._proc.is_alive()

    time.sleep(2)

    # Still spinning after two seconds -- the child did not die on us
    assert w._proc.is_alive()

    w.stop()

    # And stop() actually tore it down rather than hanging
    assert w._proc is None
    assert w._stop_queue is None
def test_WaitingAnimation__iterate(capsys):

    # _iterate runs in the child process and writes animation frames until
    # something lands in the queue. Drive it directly with an already-full
    # queue so it writes one frame and exits rather than spinning forever.
    import queue as _queue

    w = animation.WaitingAnimation(delay=0.01,num_stack=3)

    class _FullQueue:
        def empty(self):
            return False

    w._iterate(_FullQueue())

    captured = capsys.readouterr()
    assert captured.out != ""
    assert captured.out in [s for s in w._status] or captured.out.strip() != ""

def test_WaitingAnimation_start():

    w = animation.WaitingAnimation(delay=0.01)

    # Nothing running before start
    assert w._proc is None
    assert w._stop_queue is None

    w.start()
    try:
        assert w._proc is not None
        assert w._proc.is_alive()
        assert w._stop_queue is not None
    finally:
        w.stop()

def test_WaitingAnimation_stop():

    w = animation.WaitingAnimation(delay=0.01)

    # stop() before start() must be a no-op rather than an AttributeError
    w.stop()

    w.start()
    assert w._proc.is_alive()
    w.stop()

    # stop() joins the child and then closes both the process handle and the
    # queue, so the attributes are cleared. That release is what keeps repeated
    # start/stop cycles from leaking file descriptors -- we cannot query the
    # Process object afterwards precisely because it was closed.
    assert w._proc is None
    assert w._stop_queue is None

    # A second stop is still safe
    w.stop()

    # And the object can be restarted cleanly
    w.start()
    assert w._proc.is_alive()
    w.stop()
    assert w._proc is None



def test_WaitingAnimation__iterate_wraps_the_counter(capsys):
    """
    _iterate cycles through num_stack frames and then wipes the line and starts
    over. Drive it past one full cycle.
    """

    w = animation.WaitingAnimation(delay=0.001,num_stack=2)

    class _EmptyForNCalls:
        """Queue stand-in reporting empty() True for n calls, then False."""

        def __init__(self,n):
            self._n = n

        def empty(self):
            self._n -= 1
            return self._n > 0

    # More iterations than num_stack, so the counter wraps and the clear
    # string is written.
    w._iterate(_EmptyForNCalls(6))

    captured = capsys.readouterr()
    assert captured.out != ""
    assert w._clear in captured.out


def test_WaitingAnimation_stop_terminates_a_stuck_child(mocker):
    """
    If the child does not exit on its own after being told to stop, stop() has
    to terminate it -- otherwise the process (and its file descriptors) leak.
    """

    w = animation.WaitingAnimation(delay=0.01)
    w.start()

    proc = w._proc

    # Pretend the child ignored the stop signal and is still running
    mocker.patch.object(proc,"is_alive",return_value=True)
    terminate = mocker.patch.object(proc,"terminate")

    w.stop()

    terminate.assert_called_once()
    assert w._proc is None
    assert w._stop_queue is None
