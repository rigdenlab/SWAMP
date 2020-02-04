import os
import sys
import joblib
import shutil
import swamp
import argparse
import traceback
import itertools
from swamp.mr.mrjob import MrJob
from swamp.mr.mrarray import MrArray
from swamp.logger.swamplogger import SwampLogger
from swamp.command_line import check_file_exists
from swamp.library_scan.fagmentranking import FragmentRanking


def parse_arguments():
    """Parse command line arguments for mr run"""

    parser = argparse.ArgumentParser(description='SWAMP-MR: Contact assisted fragment based molecular replacement',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("id", type=str, help='Unique identifier for this MR subroutine')
    parser.add_argument("workdir", type=str, help='Working directory to perform the MR')
    parser.add_argument("mtzfile", type=check_file_exists, help='MTZ file with the reflection data')
    parser.add_argument("fastafile", type=check_file_exists, help="FASTA file with the sequence of the structure")
    parser.add_argument("conpred", type=check_file_exists, help="Residue contact prediction for the target structure")
    parser.add_argument("sspred", type=check_file_exists, help="Secondary structure prediction for the target protein")
    parser.add_argument("-nprocs", type=int, nargs="?", default=1, help="Number of parallel processors to use")
    parser.add_argument("-pdb_benchmark", type=check_file_exists, nargs="?", default=None,
                        help="PDB file with the solve structure (for benchmarking)")
    parser.add_argument("-queue_platform", type=str, nargs="?", default='sge', help="Platform to execute MR runs")
    parser.add_argument("-mtz_phases", type=check_file_exists, nargs="?", default=None,
                        help="MTZ file with phase information (for benchmarking)")
    parser.add_argument("-job_kill_time", type=int, nargs="?", default=4320, help='Maximum runtime of each MR run')
    parser.add_argument("-ncontacts_threshold", type=int, nargs="?", default=28,
                        help='Minimum no. of interhelical contacts to compute CMO of subtarget')
    parser.add_argument("-consco_threshold", type=float, nargs="?", default=0.75,
                        help='Minimum CMO between predicted and observed contacts to use search model')
    parser.add_argument('-python_interpreter', type=str, nargs='?', help='Use a given python interpreter on MR runs',
                        default=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'))
    parser.add_argument("-combine_searchmodels", action='store_true',
                        help='If set, SWAMP will try to combine search models matching different parts of the structure')
    args = parser.parse_args()

    return args


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

    This script will run the entire SWAMP MR subroutine for a given target structure. It will first scan the library
    of helical pairs and compute the CMO between the target's predicted contacts and the observed contacts for all the
    members of the library. Search models will be ranked according to their CMO, each of them will be taken to the
    MR pipeline using phaser > refmac > shelxe. If possible, SWAMP will also try to combine and place several search
    models if they are found to have a high CMO with different parts of the unkown structure.
    """

    args = parse_arguments()
    if not os.path.isdir(args.workdir):
        os.mkdir(args.workdir)

    global logger
    global centroid_dict
    global centroids

    logger = SwampLogger('SWAMP-MR')
    logger.init(logfile=os.path.join(args.workdir, "swamp_%s.debug" % args.id), use_console=True,
                console_level='info', logfile_level='debug')

    # Fragment ranking based on the centroids
    my_rank = FragmentRanking(os.path.join(args.workdir, 'ranking'), conpred=args.conpred, template_subset=centroids,
                              nthreads=args.nprocs, target_pdb_benchmark=args.pdb_benchmark, sspred=args.sspred,
                              alignment_algorithm_name='aleigen', logger=logger)
    logger.info('Using contacts to assess search model quality: matching predicted contacts with observed contacts\n')
    my_rank.rank(n_contacts_threshold=args.ncontacts_threshold)
    my_rank.rank_searchmodels(ncontacts_threshold=args.ncontacts_threshold, consco_threshold=args.consco_threshold,
                              combine_searchmodels=args.combine_searchmodels)
    shutil.rmtree(my_rank.workdir)

    # Create a new empty mr array
    my_array = MrArray(id=args.id, target_mtz=args.mtzfile, max_concurrent_nprocs=args.nprocs, target_fa=args.fastafile,
                       job_kill_time=args.job_kill_time, workdir=os.path.join(args.workdir, 'mr_array'), logger=logger,
                       platform=args.queue_platform, phased_mtz=args.mtz_phases)

    # Create MR runs for each of the ranked searchmodels
    loaded_arrangements = {}
    n_searchmodels = len(list(my_rank.ranked_searchmodels.searchmodels))
    for rank, searchmodels in enumerate(list(my_rank.ranked_searchmodels.searchmodels), 1):

        cmo = float(
            my_rank.ranked_searchmodels[my_rank.ranked_searchmodels.searchmodels == searchmodels].consco.tolist()[0])
        searchmodels = tuple([get_centroid_id(x) for x in searchmodels.split()])
        logger.debug('Loading search model no. %s/%s (CMO %s)' % (rank, n_searchmodels, cmo))

        # Only one searchmodel
        if len(searchmodels) == 1:
            combination = tuple(searchmodels)
            if not combination in loaded_arrangements.keys():
                logger.debug('Only one search model has been found in this search')
                mr_run_dir = os.path.join(workdir, 'run_1')
                mr_job = MrJob('search_%s_run_1' % rank, mr_run_dir, python_interpreter=args.python_interpreter)
                mr_job.add_searchmodel(id='1', ensemble_code=combination[0])
                loaded_arrangements[combination] = 'search_%s_run_1' % rank
                my_array.add(mr_job)
            else:
                logger.debug('%s is already used a searchmodel in the MR instance with the id '
                             '%s' % (combination, loaded_arrangements[combination]))

        # Several searchmodels must be combined
        else:
            # For each searchmodel, we get all the combinations
            logger.debug('%s search models found in this search' % len(searchmodels))
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
                        mr_job.add_searchmodel(id=idx, ensemble_code=centroid)
                    my_array.add(mr_job)

                else:
                    logger.debug('The set of search models %s is already used in the MR instance with the id '
                                 '%s' % (combination, loaded_arrangements[combination]))

    # Run task array
    my_array.run()
    my_array.append_results()
    my_array.create_result_table_outfile()
    my_array.store_pickle()


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
