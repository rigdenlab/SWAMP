"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods to manage the data in the library, and the library itself
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import os

if 'THIS_IS_READTHEDOCS' not in os.environ and "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")
