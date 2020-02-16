"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to scan the SWAMP library using contact information
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


def ScanJob(*args, **kwargs):
    """:py:obj:`~swamp.scan.scanjob.ScanJob` instance"""
    from swamp.scan.scanjob import ScanJob

    return ScanJob(*args, **kwargs)


def ScanTarget(*args, **kwargs):
    """:py:obj:`~swamp.scan.scantarget.ScanTarget` instance"""
    from swamp.scan.scantarget import ScanTarget

    return ScanTarget(*args, **kwargs)
