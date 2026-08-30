"""
A python framework for doing ancestral sequence reconstruction using
pandas dataframes and ete4 trees as the primary data structures.
"""
__author__ = "Michael J. Harms"

def _check_for_notebook():
    """
    Check whether the code is being executed in a notebook or standard
    standard python interpreter.

    Return
    ------
        string for jupyter or IPython, None for something not recognized.
    """

    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return "jupyter"   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return "IPython"   # Terminal running IPython
        else:
            return None        # Not sure what interpreter

    # Probably standard Python interpreter
    except NameError:
        return None

_in_notebook = _check_for_notebook()

# Submodules
from . import util
from . import pipeline
from . import _private
from . import quality

# Core functions for pipeline
from .quality import shrink_dataset

from .pipeline import seed_to_alignment, alignment_to_ancestors, bootstrap_reconcile

from .util import create_nicknames

from .muscle import align
from .opentree import get_df_ott, df_to_species_tree
from .ncbi.blast import recip_blast
from .raxml import find_best_model, generate_ml_tree, generate_ancestors, generate_bootstraps
from .generax import reconcile

# Input/output functions
from .io import df_from_seed
from .io import write_dataframe, read_dataframe
from .io import read_fasta_into, write_fasta, write_phy

# Topiary version
from .__version__ import __version__


# Subpackages that are expensive to import and are only needed for specific
# tasks. `draw` pulls in matplotlib, toytree, and ete4 (which itself drags in
# scipy); together these are the majority of topiary's import time and are not
# needed to run most of the pipeline. They are loaded on first attribute access
# instead, so `topiary.draw.tree(...)` still works as before.
_LAZY_SUBMODULES = ("draw",)

def __getattr__(name):
    """
    Import expensive subpackages on first access (see PEP 562).
    """

    if name in _LAZY_SUBMODULES:
        import importlib
        module = importlib.import_module(f".{name}",__name__)

        # Cache on the package so this __getattr__ is not hit again
        globals()[name] = module
        return module

    err = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(err)

def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_SUBMODULES))
