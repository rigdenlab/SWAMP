"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to prepare search models for MR
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version
import os

__version__ = version.__version__

if 'THIS_IS_READTHEDOCS' not in os.environ:

    from distutils.version import StrictVersion
    import mrbump

    if "CCP4" not in os.environ:
        raise RuntimeError("Cannot find CCP4 root directory")
    if StrictVersion(mrbump.__version__) < StrictVersion("2.0.1"):
        raise RuntimeError("MrBump must be version >= 2.0.1")


def Core(*args, **kwargs):
    """:py:obj:`~swamp.searchmodel.core.Core` instance"""
    from swamp.searchmodel.core import Core

    return Core(*args, **kwargs)


def PolyALA(*args, **kwargs):
    """:py:obj:`~swamp.searchmodel.polyala.PolyALA` instance"""
    from swamp.searchmodel.polyala import PolyALA

    return PolyALA(*args, **kwargs)


def SearchModel(*args, **kwargs):
    """:py:obj:`~swamp.searchmodel.searchmodel.SearchModel` instance"""
    from swamp.searchmodel.searchmodel import SearchModel

    return SearchModel(*args, **kwargs)
