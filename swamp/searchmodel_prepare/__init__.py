"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to rank search models according to the contact maximum overlap respect the
predicted contacts of a given target structure of interest.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import mrbump
from swamp import version
from distutils.version import StrictVersion

__version__ = version.__version__

import os

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")

if StrictVersion(mrbump.__version__) < StrictVersion("2.0.5"):
    raise RuntimeError("MrBump must be version >= 2.0.5")
