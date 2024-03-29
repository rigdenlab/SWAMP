"""SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods necessary to do molecular replacement on a given target using helical pairs.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version
import os

__version__ = version.__version__

if 'DISABLE_DEPENDENCY_CHECKS' not in os.environ:

    from distutils.version import StrictVersion

    if "CCP4" not in os.environ:
        raise RuntimeError("Cannot find CCP4 root directory")

    try:
        import dill
    except ImportError:
        raise ImportError('Dill must be installed before using SWAMP-MR')

    try:
        import prettytable
    except ImportError:
        raise ImportError('Prettytable must be installed before using SWAMP-MR')

    if "CCP4" not in os.environ:
        raise RuntimeError("Cannot find CCP4 root directory")

    if StrictVersion(prettytable.__version__) < StrictVersion("0.7.2"):
        raise RuntimeError("Prettytable must be version >= 0.7.2 to use SWAMP-MR")

    if StrictVersion(dill.__version__[:-2]) < StrictVersion("0.3"):
        raise RuntimeError("Dill must be version >= 0.3 to use SWAMP-MR")


def Mr(*args, **kwargs):
    """:py:obj:`~swamp.mr.mr.Mr` instance"""
    from swamp.mr.mr import Mr

    return Mr(*args, **kwargs)


def MrRun(*args, **kwargs):
    """:py:obj:`~swamp.mr.mrrun.MrRun` instance"""
    from swamp.mr.mrrun import MrRun

    return MrRun(*args, **kwargs)


def MrArray(*args, **kwargs):
    """:py:obj:`~swamp.mr.mrarray.MrArray` instance"""
    from swamp.mr.mrarray import MrArray

    return MrArray(*args, **kwargs)


def MrJob(*args, **kwargs):
    """:py:obj:`~swamp.mr.mrjob.MrJob` instance"""
    from swamp.mr.mrjob import MrJob

    return MrJob(*args, **kwargs)


def MrResults(*args, **kwargs):
    """:py:obj:`~swamp.mr.mrresults.MrResults` instance"""
    from swamp.mr.mrresults import MrResults

    return MrResults(*args, **kwargs)
