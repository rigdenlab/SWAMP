import os
import sys
import shutil
import joblib
import swamp
import argparse
import traceback
import pandas as pd
from swamp.utils import compress, SwampLibrary
from swamp.searchmodel import Core
from swamp.logger import SwampLogger
from swamp.command_line import check_path_exists
from swamp.clustering import Optics


def parse_arguments():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(description='SWAMP-MAKE-LIBRARY: Create an ensemble library for SWAMP',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("workdir", type=str, help='Working directory to perform ensemble clustering')
    parser.add_argument("-nprocs", type=int, nargs="?", default=1, help='Number of processors to use')
    parser.add_argument("-homologs", type=check_path_exists, nargs="?", default=None,
                        help='A file with the list of homolog structures to exclude')
    parser.add_argument("-overwrite_library", action='store_true',
                        help='If set, overwrite the SWAMP library with the new ensembles')
    parser.add_argument("-min_samples", type=int, nargs="?", default=2,
                        help='sklearn.OPTICS: no. of samples in a neighborhood for a point to be considered as a core')
    parser.add_argument("-xi", type=float, nargs="?", default=0.01,
                        help='sklearn.OPTICS: min.steepness on the reachability plot to constitute a cluster boundary')
    parser.add_argument("-cluster_method", type=str, nargs="?", default='xi',
                        help='sklearn.OPTICS: extraction method using the calculated cluster reachability')
    parser.add_argument("-max_eps", type=float, nargs="?", default=0.2,
                        help='sklearn.OPTICS: max. dist. between points to consider within neighborhood of each other')
    parser.add_argument("-eps", type=float, nargs="?", default=0.2,
                        help='sklearn.OPTICS: max. dist. between points to consider within neighborhood of each other')
    parser.add_argument("-min_cluster_size", type=int, nargs="?", default=2,
                        help='sklearn.OPTICS: min. no. of samples in a cluster')

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

    # ------------------ PARSE ARGUMENTS AND CREATE LOGGER ------------------

    args = parse_arguments()

    idx = 0
    workdir = os.path.join(args.workdir, 'SWAMP_%s' % idx)
    while os.path.isdir(workdir):
        idx += 1
        workdir = os.path.join(args.workdir, 'SWAMP_%s' % idx)
    os.mkdir(workdir)

    global logger
    logger = SwampLogger('SWAMP-MAKE-CLUSTERS')
    logger.init(logfile=os.path.join(workdir, "swamp_clustering.log"), use_console=True, console_level='info',
                logfile_level='debug')
    logger.info("Invoked with command-line:\n%s\n", " ".join(map(str, ['swamp-make-library'] + sys.argv[1:])))

    # --------------- LOAD SWAMP LIBRARY AND REMOVE HOMOLOGS ----------------

    logger.info('Loading fragment library\n')
    my_library = SwampLibrary(swamp.FRAG_PDB_DB, logger=logger)
    my_library.qscore_matrix = pd.read_csv(swamp.QSCORE_MTX_CSV, header=0, index_col=0)
    my_library.rmsd_matrix = pd.read_csv(swamp.RMSD_MTX_CSV, header=0, index_col=0)
    my_library.nalign_matrix = pd.read_csv(swamp.NALIGN_MTX_CSV, header=0, index_col=0)

    if args.homologs is not None:
        homologs_list = get_homologs(args.homologs)
        logger.info('Removing homologs from file %s => %s\n' % (args.homologs, ', '.join(homologs_list)))
        my_library.remove_homologs(homologs_list)

    # --------------------------- COMPUTE CLUSTERS --------------------------

    my_cluster = Optics(logger=logger)
    my_cluster.best_params = {'min_samples': args.min_samples, 'xi': args.xi, 'cluster_method': args.cluster_method,
                              'max_eps': args.max_eps, 'eps': args.eps, 'min_cluster_size': args.min_cluster_size}
    my_cluster.register_library(my_library)
    logger.info('Clustering fragment library to create ensembles\n')
    my_cluster.cluster()

    logger.info('Generating index of ensembles')
    my_cluster.get_cluster_dict(my_cluster.labels, inplace=True)
    my_cluster.get_centroid_dict()
    my_cluster.get_ensemble_dict(qscore_threshold=0.83, nthreads=args.nprocs)

    # ---------------------------- DUMP ENSEMBLES ---------------------------

    if args.overwrite_library:
        output_dir = swamp.ENSEMBLE_DIR
        cluster_composition_pckl = swamp.CLUSTER_COMPOSITION_PCKL
        centroid_dict_pckl = swamp.CENTROID_DICT_PCKL
    else:
        output_dir = workdir
        cluster_composition_pckl = os.path.join(workdir, 'cluster_composition.pckl')
        centroid_dict_pckl = os.path.join(workdir, 'centroid_dict.pckl')

    for ensebmle_id in my_cluster.ensemble_dict.keys():
        if my_cluster.ensemble_dict[ensebmle_id] is not None:
            # CENTROID
            centroid_fname = os.path.join(output_dir, 'centroid_%s.pdb' % ensebmle_id)

            shutil.copyfile(os.path.join(swamp.FRAG_PDB_DB, '%s.pdb' % my_cluster.centroid_dict[ensebmle_id]),
                            centroid_fname)

            # ENSEMBLE
            full_ensemble_fname = os.path.join(output_dir, 'ensemble_%s.pdb' % ensebmle_id)
            my_cluster.ensemble_dict[ensebmle_id].write_pdb(full_ensemble_fname)

            # CORE
            core_ensemble_fname = os.path.join(output_dir, 'core_%s.pdb' % ensebmle_id)
            core = Core(workdir=os.path.join(swamp.TMP_DIR, 'core_workdir'), pdbin=full_ensemble_fname, logger=logger,
                        pdbout=core_ensemble_fname)
            core.prepare()
            shutil.rmtree(core.workdir)

            # COMPRESS
            compress(full_ensemble_fname)
            os.remove(full_ensemble_fname)
            compress(centroid_fname)
            os.remove(centroid_fname)
            compress(core_ensemble_fname)
            os.remove(core_ensemble_fname)

    joblib.dump(my_cluster.composition_dict, cluster_composition_pckl, protocol=2)
    joblib.dump(my_cluster.centroid_dict, centroid_dict_pckl, protocol=2)


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
