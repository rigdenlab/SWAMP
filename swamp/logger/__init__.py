"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods to be used as a logging interface by SWAMP.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version

__version__ = version.__version__

import os

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")

_PACKAGE_PATH = os.path.join(os.environ["CCP4"], "lib", "py2", "swamp")
IDEALHELICES_DIR = os.path.join(_PACKAGE_PATH, "idealhelices")
