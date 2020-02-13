"""This is SWAMP: Solving Structures With Alpha Membrane Pairs"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version

__version__ = version.__version__

import os
import sys

TMP_DIR = os.environ["CCP4_SCR"]
PACKAGE_PATH = os.path.join(os.environ["CCP4"], "lib", "py2", "swamp")
LIBRARY = os.path.join(PACKAGE_PATH, "library", "db")
IDEALHELICES_DIR = os.path.join(LIBRARY, "idealhelices")
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

SRC_MAPALIGN = os.path.join(os.environ["CCP4"], "bin", "map_align")
SRC_ALEIGEN = os.path.join(os.environ["CCP4"], "bin", "aleigen")
SRC_WEIGENVECT = os.path.join(os.environ["CCP4"], "bin", "weigenvect")
if not os.path.isfile(SRC_ALEIGEN):
    SRC_ALEIGEN = None
if not os.path.isfile(SRC_WEIGENVECT):
    SRC_WEIGENVECT = None
if not os.path.isfile(SRC_MAPALIGN):
    SRC_MAPALIGN = None

if not os.path.isdir(LIBRARY):
    warnings.warn('SWAMP cannot find ensemble library. You can create this library using swamp-make-library')
