"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods necessary to do molecular replacement on a given target using helical pairs.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import os
from swamp import version
from distutils.version import StrictVersion

__version__ = version.__version__

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
    raise RuntimeError("Prettytable must be version >= 0.7.2")

if StrictVersion(dill.__version__[:-2]) < StrictVersion("0.3.1"):
    raise RuntimeError("Dill must be version >= 0.3.1")
