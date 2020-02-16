"""SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods to cluster the fragments present in the SWAMP library to form ensembles that
can be used as search models.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

import os
from swamp import version

__version__ = version.__version__

if 'THIS_IS_READTHEDOCS' not in os.environ:

    import scipy
    import sklearn
    from distutils.version import StrictVersion

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


def Clustering(*args, **kwargs):
    """:py:obj:`~swamp.clustering.clustering.Clustering` instance"""
    from swamp.clustering.clustering import Clustering

    return Clustering(*args, **kwargs)


def Spectral(*args, **kwargs):
    """:py:obj:`~swamp.clustering.spectral.Spectral` instance"""
    from swamp.clustering.spectral import Spectral

    return Spectral(*args, **kwargs)


def SwampAffinityPropagation(*args, **kwargs):
    """:py:obj:`~swamp.clustering.affinity_propagation.SwampAffinityPropagation` instance"""
    from swamp.clustering.swampaffinitypropagation import SwampAffinityPropagation

    return SwampAffinityPropagation(*args, **kwargs)


def Optics(*args, **kwargs):
    """:py:obj:`~swamp.clustering.optics.Optics` instance"""
    from swamp.clustering.optics import Optics

    return Optics(*args, **kwargs)
