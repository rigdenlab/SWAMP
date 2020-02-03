"""This is SWAMP: Solving Structures With Alpha Membrane Pairs"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

# ------------------ NEED TO CHANGE THIS -----------------------

SRC_MAPALIGN = "/home/filo/opt/map_align_v1/map_align/map_align"
SRC_ALEIGEN = "/home/filo/opt/AL_Eigen/aleigen"
SRC_WEIGENVECT = "/home/filo/opt/AL_Eigen/weigenvect"

# ------------------ DON'T CHANGE BELOW THIS ------------------


from swamp import version

__version__ = version.__version__

import os
import sys
import gemmi
import pyjob
import conkit
import pandas
import numpy
import Bio
import warnings

from distutils.version import StrictVersion

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")

if StrictVersion(numpy.__version__) < StrictVersion("1.16"):
    raise RuntimeError("Numpy must be version >= 1.16")

if StrictVersion(gemmi.__version__) < StrictVersion("0.3.1"):
    raise RuntimeError("Gemmi must be version >= 0.3.1")

if StrictVersion(pyjob.__version__) < StrictVersion("0.4"):
    raise RuntimeError("Pyjob must be version >= 0.4")

if StrictVersion(conkit.__version__) < StrictVersion("0.11.3"):
    raise RuntimeError("Conkit must be version >= 0.11.3")

if StrictVersion(Bio.__version__) < StrictVersion("1.74"):
    raise RuntimeError("Biopython must be version >= 1.74")

if StrictVersion(pandas.__version__) < StrictVersion("0.24"):
    raise RuntimeError("Pandas must be version >= 0.24")

TMP_DIR = os.environ["CCP4_SCR"]
PACKAGE_PATH = os.path.join(os.environ["CCP4"], "lib", "py2", "swamp")
IDEALHELICES_DIR = os.path.join(PACKAGE_PATH, "idealhelices")
LIBRARY = os.path.join(PACKAGE_PATH, "library")
FRAG_MAPALIGN_DB = os.path.join(LIBRARY, "mapalign")
FRAG_EIGEN_DB = os.path.join(LIBRARY, "eigen")
FRAG_ALEIGEN_DB = os.path.join(LIBRARY, "aleigen")
FRAG_PDB_DB = os.path.join(LIBRARY, "pdb")
ENSEMBLE_DIR = os.path.join(LIBRARY, 'ensembles')
DIST_MTX_DIR = os.path.join(LIBRARY, 'distance_mtx')
CLUSTER_COMPOSITION_PCKL = os.path.join(ENSEMBLE_DIR, 'cluster_composition.pckl')
CENTROID_DICT_PCKL = os.path.join(ENSEMBLE_DIR, 'centroid_dict.pckl')
QSCORE_MTX_CSV = os.path.join(DIST_MTX_DIR, 'qscore_mtx.csv')
RMSD_MTX_CSV = os.path.join(DIST_MTX_DIR, 'rmsd_mtx.csv')
NALIGN_MTX_CSV = os.path.join(DIST_MTX_DIR, 'nalign_mtx.csv')
if sys.hexversion >= 0x3000000:
    QSCORE_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'qscore_mtx_protocol3.pckl')
    RMSD_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'rmsd_mtx_protocol3.pckl')
    NALIGN_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'nalign_mtx_protocol3.pckl')
else:
    QSCORE_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'qscore_mtx_protocol2.pckl')
    RMSD_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'rmsd_mtx_protocol2.pckl')
    NALIGN_MTX_PCKL = os.path.join(DIST_MTX_DIR, 'nalign_mtx_protocol2.pckl')
if not os.path.isfile(SRC_ALEIGEN):
    SRC_ALEIGEN = None
if not os.path.isfile(SRC_MAPALIGN):
    SRC_MAPALIGN = None

if not os.path.isdir(LIBRARY):
    warnings.warn('SWAMP cannot find ensemble library. You can create this library using swamp-make-library')
