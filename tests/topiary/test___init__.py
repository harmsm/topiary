import pytest

from topiary.__init__ import _check_for_notebook


def test__check_for_notebook(mocker):

    # _check_for_notebook looks up get_ipython in the globals of the module it
    # was defined in. Note that `topiary` and `topiary.__init__` are two
    # distinct module objects, so patch the function's own globals rather than
    # guessing which alias holds them.
    func_globals = _check_for_notebook.__globals__

    def _shell_named(name):
        """Return a fake get_ipython whose result's class is called `name`."""
        fake_shell = type(name,(object,),{})()
        return lambda: fake_shell

    # No get_ipython at all -> NameError branch -> None. This is the normal
    # case when topiary is imported from a plain interpreter.
    assert "get_ipython" not in func_globals
    assert _check_for_notebook() is None

    mocker.patch.dict(func_globals,
                      {"get_ipython":_shell_named("ZMQInteractiveShell")})
    assert _check_for_notebook() == "jupyter"

    mocker.patch.dict(func_globals,
                      {"get_ipython":_shell_named("TerminalInteractiveShell")})
    assert _check_for_notebook() == "IPython"

    # Some shell we do not recognize
    mocker.patch.dict(func_globals,
                      {"get_ipython":_shell_named("SomethingElseEntirely")})
    assert _check_for_notebook() is None
