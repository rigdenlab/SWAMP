"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods related to molecular replacement and search model preparation
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import os
from swamp import version

__version__ = version.__version__

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")
