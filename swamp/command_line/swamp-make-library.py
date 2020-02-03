import os
import sys
import shutil
import joblib
import swamp
import argparse
import traceback
import pandas as pd
from swamp.command_line import check_file_exists
from swamp.clustering.optics import Optics
from swamp.logger.swamplogger import SwampLogger
from swamp.library_tools.swamplibrary import SwampLibrary
from swamp.searchmodel_prepare.core import Core


def parse_arguments():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(description='SWAMP-MAKE-LIBRARY: Create an ensemble library for SWAMP',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("workdir", type=str, help='Working directory to perform ensemble clustering')
    parser.add_argument("-nprocs", type=int, nargs="?", default=1, help='Number of processors to use')
    parser.add_argument("-homologs", type=check_file_exists, nargs="?", default=None,
                        help='A file with the list of homolog structures to exclude')
    parser.add_argument("-overwrite_library", action='store_true',
                        help='If set, overwrite the SWAMP library with the new ensembles')
    args = parser.parse_args()

    return args


def get_homologs(homologs_fname):
    """Extract the list of homolog structures from a list file

    :param homologs_fname: file name with the list of homologs
    :type homologs_fname: str
    :returns homologs_list containing the pdb codes of the homolog structures in the input file
    :rtype tuple
    """

    homologs_list = []
    with open(homologs_fname, 'r') as fhandle:
        for line in fhandle:
            homologs_list.append(line.rstrip()[:-2].lower())
    return tuple(homologs_list)


def main():
    """Execute a swamp make clusters

    This script will create a SWAMP library based on the user's input. It can be used to update the library with new
    ensembles, to create a library excluding certain structures (homolog removal for benchmarking) or to create
    ensembles using a different clustering algorithm
    """

    args = parse_arguments()
    if not os.path.isdir(args.workdir):
        os.mkdir(args.workdir)

    global logger
    logger = SwampLogger('SWAMP-MAKE-CLUSTERS')
    logger.init(logfile=os.path.join(args.workdir, "swamp_clustering.log"), use_console=True, console_level='info',
                logfile_level='debug')

    # Load library
    logger.info('Loading fragment library\n')
    my_library = SwampLibrary(swamp.FRAG_PDB_DB, logger=logger)
    my_library.qscore_matrix = pd.read_csv(swamp.QSCORE_MTX_CSV, header=0, index_col=0)
    my_library.rmsd_matrix = pd.read_csv(swamp.RMSD_MTX_CSV, header=0, index_col=0)
    my_library.nalign_matrix = pd.read_csv(swamp.NALIGN_MTX_CSV, header=0, index_col=0)

    # Load homologs and remove them from the library
    if args.homologs is not None:
        homologs_list = get_homologs(args.homologs)
        logger.info('Removing homologs from file %s => %s\n' % (args.homologs, ', '.join(homologs_list)))
        my_library.remove_homologs(homologs_list)

    # Make clusters
    my_cluster = Optics(logger=logger)
    my_cluster.best_params = {'min_samples': 2, 'xi': 0.01, 'cluster_method': 'xi', 'max_eps': 0.2, 'eps': 0.2,
                              'min_cluster_size': 2}
    my_cluster.register_library(my_library)
    logger.info('Clustering fragment library to create ensembles\n')
    my_cluster.cluster()

    # Process data from clusters
    logger.info('Generating index of ensembles')
    my_cluster.get_cluster_dict(my_cluster.labels, inplace=True)
    my_cluster.get_centroid_dict()
    my_cluster.get_ensemble_dict(qscore_threshold=0.83, nthreads=35)

    if args.overwrite_library:

        # Dump the ensembles into pdb files
        for ensebmle_id in my_cluster.ensemble_dict.keys():
            if my_cluster.ensemble_dict[ensebmle_id] is not None:
                full_ensemble_fname = os.path.join(swamp.TMP_DIR, 'ensemble_%s.pdb' % ensebmle_id)
                core_ensemble_fname = os.path.join(swamp.ENSEMBLE_DIR, 'ensemble_%s.pdb' % ensebmle_id)
                my_cluster.ensemble_dict[ensebmle_id].write_pdb(full_ensemble_fname)
                core = Core(workdir=os.path.join(swamp.TMP_DIR, 'core_workdir'), pdbin=full_ensemble_fname,
                            pdbout=core_ensemble_fname, logger=logger)
                core.prepare()
                shutil.rmtree(core.workdir)

        # Save the cluster information
        joblib.dump(my_cluster.composition_dict, swamp.CLUSTER_COMPOSITION_PCKL, protocol=2)
        joblib.dump(my_cluster.centroid_dict, swamp.CENTROID_DICT_PCKL, protocol=2)


if __name__ == "__main__":
    logger = None

    try:
        main()
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)