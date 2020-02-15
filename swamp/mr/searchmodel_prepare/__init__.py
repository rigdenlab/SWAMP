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
