"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to search the SWAMP library using contact information. The contact maximum \
overlap (CMO) will be calculated between the predicted contacts and the observed contacts.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version

__version__ = version.__version__

import os

if 'THIS_IS_READTHEDOCS' not in os.environ:

    if "CCP4" not in os.environ:
        raise RuntimeError("Cannot find CCP4 root directory")

    try:
        import joblib
    except ImportError:
        raise ImportError('Joblib must be installed before using SWAMP-SCAN')


def SearchJob(*args, **kwargs):
    """:py:obj:`~swamp.search.searchjob.SearchJob` instance"""
    from swamp.search.searchjob import SearchJob

    return SearchJob(*args, **kwargs)


def SearchTarget(*args, **kwargs):
    """:py:obj:`~swamp.search.searchtarget.SearchTarget` instance"""
    from swamp.search.searchtarget import SearchTarget

    return SearchTarget(*args, **kwargs)
