import os
import sys
import joblib
import swamp
import argparse
import traceback
import itertools
from swamp.logger import SwampLogger
from swamp.search import SearchTarget
from swamp.mr import MrJob, MrArray, MrResults
from swamp.command_line import check_path_exists


def create_argument_parser():
    """Create a parser for the command line arguments used in swamp-mr"""

    parser = argparse.ArgumentParser(description='SWAMP-MR: Contact assisted fragment based molecular replacement',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("id", type=str, help='Unique identifier for this MR subroutine')
    parser.add_argument("mtzfile", type=check_path_exists, help='MTZ file with the reflection data')
    parser.add_argument("fastafile", type=check_path_exists, help="FASTA file with the sequence of the structure")
    parser.add_argument("conpred", type=check_path_exists, help="Residue contact prediction for the target structure")
    parser.add_argument("sspred", type=check_path_exists, help="Secondary structure prediction for the target protein")
    parser.add_argument("-workdir", type=str, nargs='?', default=None, help='Working directory to perform MR')
    parser.add_argument("-max_concurrent_procs", type=int, nargs="?", default=1, help="Max no. of concurrent processes")
    parser.add_argument("-max_array_size", type=int, nargs="?", default=None,
                        help='Maximum allowed array size to be submitted to the job scheduler')
    parser.add_argument("-pdb_benchmark", type=check_path_exists, nargs="?", default=None,
                        help="PDB file with the solve structure (for benchmarking)")
    parser.add_argument("-platform", type=str, nargs="?", default='sge',
                        help="Platform to execute MR runs (default sge)")
    parser.add_argument("-environment", type=str, nargs="?", default=None,
                        help="Select a environment for execution of HPC tasks (default None)")
    parser.add_argument("-queue", type=str, nargs="?", default=None,
                        help="Queue name to send jobs in the HPC (default None)")
    parser.add_argument("-mtz_phases", type=check_path_exists, nargs="?", default=None,
                        help="MTZ file with phase information (for benchmarking)")
    parser.add_argument("-job_kill_time", type=int, nargs="?", default=4320, help='Maximum runtime of each MR run')
    parser.add_argument("-ncontacts_threshold", type=int, nargs="?", default=28,
                        help='Minimum no. of interhelical contacts to compute CMO of subtarget')
    parser.add_argument("-consco_threshold", type=float, nargs="?", default=0.75,
                        help='Minimum CMO between predicted and observed contacts to use search model')
    parser.add_argument('-python_interpreter', type=str, nargs='?', help='Indicate a python interpreter for MR runs',
                        default=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'))
    parser.add_argument("-combine_searchmodels", action='store_true',
                        help='If set, combine search models matching different parts of the structure')
    parser.add_argument("-use_centroids", action='store_true',
                        help='Centroids used as search models as well (only if not combining search models)')
    parser.add_argument("-use_cores", action='store_true',
                        help='Core ensembles used as search models as well (only if not combining search models)')
    return parser


def get_centroid_id(frag_id):
    """Get the cluster id for a given centroid

    :param frag_id: identifier of the fragment of interest
    :type frag_id: int
    :returns cluster id where the fragment is found
    :rtype str
    """

    global centroid_dict
    for clst_id in centroid_dict.keys():
        if frag_id == centroid_dict[clst_id]:
            return clst_id


def main():
    """Execute a swamp mr run

    This script will run the entire SWAMP MR subroutine for a given target structure. It will first search the library
    of helical pairs and compute the CMO between the target's predicted contacts and the observed contacts for all the
    members of the library. Search models will be ranked according to their CMO, each of them will be taken to the
    MR pipeline using phaser > refmac > shelxe. If possible, SWAMP will also try to combine and place several search
    models if they are found to have a high CMO with different parts of the unkown structure.
    """

    # ------------------ PARSE ARGUMENTS AND CREATE LOGGER ------------------

    parser = create_argument_parser()
    args = parser.parse_args()

    if args.workdir is not None:
        if os.path.isdir(args.workdir):
            raise FileExistsError('Working directory already exists!')
        swamp_workdir = args.workdir
    else:
        idx = 0
        swamp_workdir = os.path.join(os.getcwd(), 'SWAMP_%s' % idx)
        while os.path.isdir(swamp_workdir):
            idx += 1
            swamp_workdir = os.path.join(os.getcwd(), 'SWAMP_%s' % idx)

    os.mkdir(swamp_workdir)
    logfile = os.path.join(swamp_workdir, "swamp_%s.debug" % args.id)
    swamp_search_dir = os.path.join(swamp_workdir, 'swamp_search')
    swamp_mr_dir = os.path.join(swamp_workdir, 'swamp_mr')

    global logger
    global centroid_dict
    global centroids

    logger = SwampLogger('SWAMP-MR')
    logger.init(logfile=logfile, use_console=True,
                console_level='info', logfile_level='debug')
    logger.info("Invoked with command-line:\n%s\n", " ".join(map(str, ['swamp-mr'] + sys.argv[1:])))

    # ------------------ SCAN LIBRARY OF SEARCH MODELS USING CONTACTS ------------------

    my_rank = SearchTarget(swamp_search_dir, conpred=args.conpred, template_subset=centroids, sspred=args.sspred,
                           alignment_algorithm_name='aleigen', n_contacts_threshold=args.ncontacts_threshold,
                           nthreads=args.max_concurrent_procs, target_pdb_benchmark=args.pdb_benchmark,
                           logger=logger, platform=args.platform, python_interpreter=args.python_interpreter,
                           queue_environment=args.environment, queue_name=args.queue)
    logger.info('Using contacts to assess search model quality: matching predicted contacts with observed contacts\n')
    my_rank.search()
    my_rank.rank(consco_threshold=args.consco_threshold, combine_searchmodels=args.combine_searchmodels)
    if args.pdb_benchmark is not None:
        joblib.dump(my_rank.results, os.path.join(swamp_search_dir, 'swamp-search.pckl'))

    # ------------------ CREATE MR TASK ARRAY AND LOAD INDIVIDUAL MR JOBS ------------------

    my_array = MrArray(id=args.id, target_mtz=args.mtzfile, max_concurrent_nprocs=args.max_concurrent_procs,
                       target_fa=args.fastafile, queue_environment=args.environment, max_array_size=args.max_array_size,
                       job_kill_time=args.job_kill_time, workdir=swamp_mr_dir, logger=logger, queue_name=args.queue,
                       platform=args.platform, phased_mtz=args.mtz_phases)

    loaded_arrangements = {}
    n_searchmodels = len(list(my_rank.ranked_searchmodels.searchmodels))
    for rank, searchmodels in enumerate(list(my_rank.ranked_searchmodels.searchmodels), 1):

        cmo = float(
            my_rank.ranked_searchmodels[my_rank.ranked_searchmodels.searchmodels == searchmodels].consco.tolist()[0])
        searchmodels = tuple([get_centroid_id(x) for x in searchmodels.split()])
        logger.debug('Loading search model no. %s/%s (CMO %s)' % (rank, n_searchmodels, cmo))
        workdir = os.path.join(my_array.workdir, 'search_%s' % rank)

        # ------------------ ONLY ONE SEARCH MODEL ------------------

        if len(searchmodels) == 1:
            combination = tuple(searchmodels)
            if combination not in loaded_arrangements.keys():

                logger.debug('Only one search model has been found in this search: %s' % combination)
                mr_run_dir = os.path.join(workdir, 'run_1')
                mr_job = MrJob('search_%s_run_1' % rank, mr_run_dir, python_interpreter=args.python_interpreter)
                mr_job.add_searchmodel(id='1', ensemble_code=combination[0], mod='polyala')
                loaded_arrangements[combination] = 'search_%s_run_1' % rank
                my_array.add(mr_job)

                if args.use_centroids:
                    logger.debug('Creating MR job for centroid at search_%s_run_2' % rank)
                    mr_run_dir = os.path.join(workdir, 'run_2')
                    mr_job = MrJob('search_%s_run_2' % rank, mr_run_dir, python_interpreter=args.python_interpreter)
                    mr_job.add_searchmodel(id='1', ensemble_code=combination[0], model='centroid', mod='polyala')
                    my_array.add(mr_job)

                if args.use_cores:
                    logger.debug('Creating MR job for core ensemble at search_%s_run_3' % rank)
                    mr_run_dir = os.path.join(workdir, 'run_3')
                    mr_job = MrJob('search_%s_run_3' % rank, mr_run_dir, python_interpreter=args.python_interpreter)
                    mr_job.add_searchmodel(id='1', ensemble_code=combination[0], model='core', mod='polyala')
                    my_array.add(mr_job)

            else:
                logger.debug('%s is already used a searchmodel in the MR instance with the id '
                             '%s' % (combination, loaded_arrangements[combination]))

        # ------------------ COMBINE MULTIPLE SEARCH MODELS ------------------

        else:
            # For each searchmodel, we get all the combinations
            logger.debug('%s search models found in this search: %s' % (len(searchmodels), searchmodels))
            search_list = []
            for nsearch in range(1, len(searchmodels) + 1):
                search_list += list(itertools.combinations(searchmodels, nsearch))
            logger.debug('Loading %s possible arrangements' % len(search_list))

            # For each search model combination we create separate runs
            for nrun, combination in enumerate(search_list, 1):

                if combination not in loaded_arrangements.keys():
                    job_id = 'search_%s_run_%s' % (rank, nrun)
                    mr_run_dir = os.path.join(workdir, 'run_%s' % nrun)
                    loaded_arrangements[combination] = mr_run_dir
                    logger.debug('Processing arrangement %s/%s (%s)' % (nrun, len(search_list), mr_run_dir))
                    mr_job = MrJob(job_id, mr_run_dir, python_interpreter=args.python_interpreter)
                    for idx, centroid in enumerate(combination, 1):
                        mr_job.add_searchmodel(id=idx, ensemble_code=centroid, mod='polyala')
                    my_array.add(mr_job)

                else:
                    logger.debug('The set of search models %s is already used in the MR instance with the id '
                                 '%s' % (combination, loaded_arrangements[combination]))

    # ------------------ SUBMIT MR TASK ARRAY FOR EXECUTION ------------------

    my_array.run()
    results = MrResults(os.path.join(swamp_workdir), logger=logger)
    results.report_results()


if __name__ == "__main__":

    logger = None
    centroid_dict = joblib.load(swamp.CENTROID_DICT_PCKL)
    centroids = list(centroid_dict.values())

    try:
        main()
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
