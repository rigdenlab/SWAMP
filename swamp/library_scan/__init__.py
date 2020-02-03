"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to prepare search models before molecular replacement takes place.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version

__version__ = version.__version__

import os

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")
