"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods to cluster the fragments present in the SWAMP library to form ensembles that
can be used as search models.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import scipy
import sklearn
from swamp import version
from distutils.version import StrictVersion

__version__ = version.__version__

import os

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")

if StrictVersion(sklearn.__version__) < StrictVersion("0.21"):
    raise RuntimeError("Sklearn must be version >= 0.21")

if StrictVersion(scipy.__version__) < StrictVersion("1.3.1"):
    raise RuntimeError("Scipy must be version >= 1.3.1")

try:
    import statistics
except ImportError as e:
    raise ImportError(e)