import os
import swamp
import joblib
import conkit.io
import itertools
import pandas as pd
from pyjob import TaskFactory
from swamp.scan import ScanJob
from swamp.logger import SwampLogger
from swamp.utils import TargetSplit, SwampLibrary, renumber_hierarchy


class ScanTarget(object):
    """Class to scan the library and rank search models according to their CMO with a given target.

    Using CMO alignment tools determine the best fragments in the library to be used as search models. The target
    will be split into several subtargets (one for each helical pair with enough interhelical contact information,
    and several scan jobs will take place (one for each subtarget).

    :param str workdir: working directory where the scan tasks will be executed
    :param str conpred: contact prediction file of the target
    :param str sspred: secondary structure prediction file of the target (must be topcons file)
    :param str conformat: format of the contact prediction file for the target
    :param str nthreads: number of parallel threads to use for CMO calculations (default: 1)
    :param list template_subset: set of templates to be used rather than the full fragment library (deafult: None)
    :param str target_pdb_benchmark: target's pdb file for benchmark purposes (default: None)
    :param str alignment_algorithm_name: algorithm used for CMO calculation (default: 'mapalign')
    :param `swamp.logger.swamplogger.SwampLogger` logger: logging interface for the scan (default: None)
    :param int n_contacts_threshold: min. no. of interhelical contacts to include a subtarget in the scan (default: 28)
    :param str platform: queueing system used in the HPC where the array will be executed (default 'sge')
    :param str queue_name: name of the HPC qeue where the tasks should be sent (default None)
    :param str queue_environment: name of the HPC queue environment where the tasks should be sent (default None)
    :ivar str shell_interpreter: location of the shell interpreter to be used for task execution (default '/bin/bash')
    :ivar bool error: True if errors have occurred at some point on the pipeline

    :example:
    >>> from swamp.library.scan.scantarget import ScanTarget
    >>> my_rank = ScanTarget('<workdir>', '<conpred>', '<sspred>')
    >>> my_rank.scan()
    >>> my_rank.rank()
    """

    def __init__(self, workdir, conpred, sspred, conformat="psicov", nthreads=1, template_subset=None, logger=None,
                 target_pdb_benchmark=None, alignment_algorithm_name='mapalign', n_contacts_threshold=28,
                 python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'),
                 platform='sge', queue_name=None, queue_environment=None):
        self._workdir = workdir
        self._conpred = conpred
        self._sspred = sspred
        self._conformat = conformat
        self._nthreads = nthreads
        self._target_pdb_benchmark = target_pdb_benchmark
        self._alignment_algorithm_name = alignment_algorithm_name
        self._template_subset = template_subset
        self._n_contacts_threshold = n_contacts_threshold
        self._platform = platform
        self._queue_name = queue_name
        self._queue_environment = queue_environment
        self._shell_interpreter = '/bin/bash'
        self._python_interpreter = python_interpreter
        self._con_precision_dict = None
        self._scan_pickle_dict = None
        self._scripts = None
        self._results = None
        self._make_workdir()

        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(logfile=os.path.join(self.workdir, "swamp_rank.log"), use_console=True,
                             console_level='info',
                             logfile_level='debug')
        else:
            self._logger = logger
        self._target = TargetSplit(workdir=self.workdir, conpred=self.conpred, sspred=self.sspred, logger=self.logger,
                                   conformat=self.conformat, pdb_benchmark=self.target_pdb_benchmark)

        if self.target_pdb_benchmark is not None and not os.path.isdir(swamp.FRAG_PDB_DB):
            self.logger.warning('PDB benchmark was requested but %s PDB library was not found!' % swamp.FRAG_PDB_DB)
            self.target_pdb_benchmark = None

    def __repr__(self):

        return '{}(workdir="{_workdir}", conpred="{_conpred}", sspred="{_sspred}", conformat="{_conformat}", ' \
               'nthreads="{_nthreads}", template_subset="{_template_subset}", logger="{_logger}", ' \
               'target_pdb_benchmark="{_target_pdb_benchmark}", alignment_algorithm_name="{_alignment_algorithm_name}' \
               '")'.format(self.__class__.__name__, **self.__dict__)

    # ------------------ Properties ------------------

    @property
    def scan_header(self):
        """Abstract property to store the wrapper header for the logger"""

        return """**********************************************************************
********************            SWAMP SCAN            ****************
**********************************************************************

"""

    @property
    def ranked_searchmodels(self):
        """Property to store ranked combinations of searchmodels"""
        return self._ranked_searchmodels

    @ranked_searchmodels.setter
    def ranked_searchmodels(self, value):
        self._ranked_searchmodels = value

    @property
    def n_contacts_threshold(self):
        return self._n_contacts_threshold

    @n_contacts_threshold.setter
    def n_contacts_threshold(self, value):
        self._n_contacts_threshold = value

    @property
    def queue_environment(self):
        return self._queue_environment

    @queue_environment.setter
    def queue_environment(self, value):
        self._queue_environment = value

    @property
    def queue_name(self):
        return self._queue_name

    @queue_name.setter
    def queue_name(self, value):
        self._queue_name = value

    @property
    def shell_interpreter(self):
        return self._shell_interpreter

    @shell_interpreter.setter
    def shell_interpreter(self, value):
        self._shell_interpreter = value

    @property
    def python_interpreter(self):
        return self._python_interpreter

    @python_interpreter.setter
    def python_interpreter(self, value):
        self._python_interpreter = value

    @property
    def platform(self):
        return self._platform

    @platform.setter
    def platform(self, value):
        self._platform = value

    @property
    def scripts(self):
        return self._scripts

    @scripts.setter
    def scripts(self, value):
        self._scripts = value

    @property
    def logger(self):
        """Property to store the logger interface :obj:`~swamp.logger.swamplogger.SwampLogger`"""
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def con_precision_dict(self):
        """Property to store the precision of the target's predicted contacts (only used when pdb_benchmark is set)"""
        return self._con_precision_dict

    @con_precision_dict.setter
    def con_precision_dict(self, value):
        self._con_precision_dict = value

    @property
    def target(self):
        """Target splitter :obj:`~swamp.scan.targetsplit.SplitTarget`"""
        return self._target

    @target.setter
    def target(self, value):
        self._target = value

    @property
    def template_subset(self):
        return self._template_subset

    @template_subset.setter
    def template_subset(self, value):
        self._template_subset = value

    @property
    def alignment_algorithm_name(self):
        return self._alignment_algorithm_name

    @alignment_algorithm_name.setter
    def alignment_algorithm_name(self, value):
        self._alignment_algorithm_name = value

    @property
    def nthreads(self):
        return self._nthreads

    @nthreads.setter
    def nthreads(self, value):
        self._nthreads = value

    @property
    def conformat(self):
        return self._conformat

    @conformat.setter
    def conformat(self, value):
        self._conformat = value

    @property
    def sspred(self):
        return self._sspred

    @sspred.setter
    def sspred(self, value):
        self._sspred = value

    @property
    def conpred(self):
        return self._conpred

    @conpred.setter
    def conpred(self, value):
        self._conpred = value

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    @property
    def target_pdb_benchmark(self):
        return self._target_pdb_benchmark

    @target_pdb_benchmark.setter
    def target_pdb_benchmark(self, value):
        self._target_pdb_benchmark = value

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def scan_pickle_dict(self):
        return self._scan_pickle_dict

    @scan_pickle_dict.setter
    def scan_pickle_dict(self, value):
        self._scan_pickle_dict = value

    @property
    def _tmp_cmap(self):
        return os.path.join(self.workdir, "tmp_cmap_{}.map")

    @property
    def _scan_workdir(self):
        return os.path.join(self.workdir, "scan_{}")

    @property
    def _tmp_pdb(self):
        if self.target_pdb_benchmark is not None:
            return os.path.join(self.workdir, 'tmp_strct_{}.pdb')
        else:
            return None

    @property
    def template_library(self):
        """Location of the template library to be used during the CMO scan"""

        if self.alignment_algorithm_name == 'aleigen':
            return swamp.FRAG_ALEIGEN_DB
        elif self.alignment_algorithm_name == 'mapalign':
            return swamp.FRAG_MAPALIGN_DB
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    @property
    def library_format(self):
        """Dictionary specifyinh the template library format"""

        if self.alignment_algorithm_name == 'aleigen':
            return 'aleigen'
        elif self.alignment_algorithm_name == 'mapalign':
            return 'mapalign'
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    @property
    def _other_task_info(self):
        """Dictionary with extra arguments for :py:obj:pyjob.TaskFactory"""

        info = {'directory': self.workdir, 'shell': self.shell_interpreter}

        if self.platform == 'local':
            info['processes'] = self.nthreads
        else:
            info['max_array_size'] = self.nthreads
        if self.queue_environment is not None:
            info['environment'] = self.queue_environment
        if self.queue_name is not None:
            info['queue'] = self.queue_name

        return info

    @property
    def _column_reference(self):
        """A list of fields of the result table for each of the CMO algorithms"""

        if self.alignment_algorithm_name == 'aleigen':
            return ["SUBTRGT_RANK", "SUBTRGT_ID", "N_CON_MAP_A", "MAP_A", "MAP_B", "CON_SCO", "C1", "C2", "CMO",
                    "ALI_LEN", "QSCORE", "RMSD", "SEQ_ID", "N_ALIGN"]
        elif self.alignment_algorithm_name == 'mapalign':
            return ["SUBTRGT_RANK", "SUBTRGT_ID", "N_CON_MAP_A", "MAP_A", "MAP_B", "CON_SCO", "GAP_SCO",
                    "TOTAL_SCO", "ALI_LEN", "QSCORE", "RMSD", "SEQ_ID", "N_ALIGN"]
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    # ------------------ Hidden methods ------------------

    def _make_workdir(self):
        """Create the working directory"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    def _create_scripts(self):
        """Create a list with all the :obj:`pyjob.Script` instances that need to be executed. There is one instance
        for each of the subtargets that pass the no. interhelical contact threshoold"""

        self.scripts = []
        self.con_precision_dict = {}
        self.scan_pickle_dict = {}

        for idx, subtarget in enumerate(self.target.ranked_subtargets, 1):
            self.logger.info("%s interhelical contacts found for subtarget %s" % (subtarget.ncontacts, idx))

            if subtarget.ncontacts >= self.n_contacts_threshold:

                self.logger.info('The no. of contacts passes the threshold. Creating a task array for subtarget.')
                conkit.io.write(fname=self._tmp_cmap.format(idx), format=self.library_format, hierarchy=subtarget)

                if self.target_pdb_benchmark is not None:
                    renumber_hierarchy(self.target.subtargets_pdb[subtarget.id])
                    self.target.subtargets_pdb[subtarget.id].write_pdb(self._tmp_pdb.format(idx))
                    perfect_contacts = conkit.io.read(self._tmp_pdb.format(idx), 'pdb').top_map
                    subtarget.sequence = perfect_contacts.sequence.deepcopy()
                    subtarget.set_sequence_register()
                    precision = subtarget.match(perfect_contacts).precision
                    self.con_precision_dict[subtarget.id] = precision

                scanner = ScanJob(**self._scan_info(idx))

                self.scripts.append(scanner.script)
                self.scan_pickle_dict[scanner.pickle_fname] = subtarget

            else:
                self.logger.info("No. of contacts below the %s threshold. "
                                 "Subtarget will not be used in the scan." % self.n_contacts_threshold)

        self.logger.info('%s subtargets will be used in the scan. Creating task now.' % len(self.scan_pickle_dict))

    def _scan_info(self, idx):
        """Create a dictionary with the arguments for the library scan

        :param idx: the index of the scan job
        :type idx: int
        :returns dictionary wuth the arguments to use for the library scan
        :rtype dict
        """

        info = {'workdir': self._scan_workdir.format(idx), 'pdb_library': swamp.FRAG_PDB_DB,
                'query': self._tmp_cmap.format(idx), 'algorithm': self.alignment_algorithm_name,
                'template_subset': self.template_subset, 'python_interpreter': self.python_interpreter,
                'template_library': self.template_library, 'library_format': self.library_format,
                'logger': self.logger, 'con_format': self.library_format, 'id': idx}

        if self.target_pdb_benchmark is not None:
            info['query_pdb_benchmark'] = self._tmp_pdb.format(idx),

        return info

    def _make_dataframe(self, results, **kwargs):
        """Convert the results list into a :obj:`pd.Dataframe`

        :param results: Nested list with the results of each of the CMO scans
        :type results: list
        :param kwargs: Passed to :obj:`~swamp.scan.fragrank._get_best_alignment`
        :return: no value
        :rtype: None
        """

        self.results = pd.DataFrame(results)
        self.results.columns = self._column_reference
        self.results["CENTROID_ID"] = [os.path.basename(fname).split('.')[0] for fname in self.results.MAP_B.to_list()]
        self._get_best_alignment(**kwargs)

    def _get_best_alignment(self, select_by=("CENTROID_ID", "SUBTRGT_ID")):
        """Fill the results with the optimal CMO alignment between helical pairs. This methof considers inverted
         fragments and takes highest CMO score)

         :param select_by: indicate the field by which the optimal alignment will be determined
         :type select_by: list, tuple
         :returns nothing
         :rtype None
         """

        if self.results is None:
            return

        else:
            self.results.CENTROID_ID = self.results.CENTROID_ID.apply(lambda x: SwampLibrary._get_unique_frag_id(x))
            self.results["max_score"] = self.results.groupby(list(select_by), sort=False)["CON_SCO"].transform(max)
            self.results = self.results[self.results.CON_SCO == self.results.max_score]
            self.results.drop("max_score", 1, inplace=True)
            self.results.drop_duplicates(subset=list(select_by), inplace=True)

    # ------------------ Some general methods ------------------

    def recover_results(self):
        """Recover the results from all the :obj:`swamp.library.scan.scanjob` in the :obj:`pyjob.TaskFactory`

        :returns results: a list with the results obtained for the contact map alignment across all subtargets
        :rtype list
        """

        results = []
        for pickle_fname, subtarget in zip(self.scan_pickle_dict.keys(), self.scan_pickle_dict.values()):

            if os.path.isfile(pickle_fname):
                self.logger.debug('Retrieving results from %s' % pickle_fname)
                current_job_results = joblib.load(pickle_fname)
                for result in current_job_results:
                    results.append([self.target.ranked_subtargets.index(subtarget) + 1, subtarget.id,
                                    subtarget.ncontacts] + result)

            else:
                self.logger.warning('Cannot find pickle file! %s' % pickle_fname)

        return results

    def scan(self):
        """Scan the library and compute the CMO between the observed contacts and the predicted contacts of the target.

        This method will run a :obj:`` instance that contains one scan job for each of the subtargets that met the
        no. contacts threshold.
        """

        self.logger.info(self.scan_header)
        self.logger.info("Splitting the target into sets of contacting helical pairs")
        self.target.split()
        if self.target.error:
            self.logger.warning('Previous errors prevent scanning the library with the target contacts!')
            return

        self.logger.info('Creating a list of jobs to scan the library using contacts.')
        self._create_scripts()
        self.logger.info('Sending jobs now.')

        with TaskFactory(self.platform, tuple(self.scripts), **self._other_task_info) as task:
            task.name = 'swamp_scan'
            task.run()
            self.logger.info('Waiting for workers...')

        self.logger.info('All scan tasks have been completed! Retrieving results')
        results = self.recover_results()
        self._make_dataframe(results)

    def rank(self, consco_threshold=0.75, combine_searchmodels=False):
        """Get fragments groupings ranked by their combined scores. Takes in consideration same combination of fragments
         may appear more than once with fragments in different subtargets (then the combined_consco is the maximum
         possible)

         :param consco_threshold: contact score threshold to consider a CMO alignment valid (default 0.75)
         :type consco_threshold: float
         :param combine_searchmodels: if True combine search models matching different subtargets (default False)
         :type combine_searchmodels: bool
         """

        if self.results is None:
            self.logger.error('Need to rank the fragments first!')
            return

        # Get the valid searchmodels (subtarget with enough contacts and fragment with enough consco)
        valid_searchmodels = self.results.loc[self.results.CON_SCO >= consco_threshold,]
        if len(set(list(valid_searchmodels.SUBTRGT_ID))) == 0:
            self.logger.warning('None of the search models meet CMO criteria!')
            return
        elif not combine_searchmodels:
            self.logger.info('%s search models meet CMO criteria' % len(set(valid_searchmodels.CENTROID_ID.tolist())))
            self.ranked_searchmodels = pd.DataFrame()
            self.ranked_searchmodels['searchmodels'] = valid_searchmodels.CENTROID_ID
            self.ranked_searchmodels['consco'] = valid_searchmodels.CON_SCO
            self.ranked_searchmodels.sort_values(by='consco', inplace=True, ascending=False)
            return
        elif len(set(list(valid_searchmodels.SUBTRGT_ID))) == 1:
            self.logger.warning('Only subtarget %s has valid searchmodels!' % list(valid_searchmodels.SUBTRGT_ID)[0])
            self.ranked_searchmodels = pd.DataFrame()
            self.ranked_searchmodels['searchmodels'] = valid_searchmodels.CENTROID_ID
            self.ranked_searchmodels['consco'] = valid_searchmodels.CON_SCO
            self.ranked_searchmodels.sort_values(by='consco', inplace=True, ascending=False)
            return

        self.logger.info('Retrieving ranked combinations of search models...')

        # Get a dictionary with the valid fragments for each subtarget
        all_valid_searchmodels = {}
        for subtarget in set(valid_searchmodels.SUBTRGT_ID):
            all_valid_searchmodels[subtarget] = list(
                valid_searchmodels[valid_searchmodels.SUBTRGT_ID == subtarget].CENTROID_ID)

        # All valid combinations between valid searchmodels
        all_valid_combinations = list(itertools.product(*list(all_valid_searchmodels.values())))

        # Make a dictionary with the searchmodel combination and the combined consco
        searchmodel_combinations = {}
        for combination in all_valid_combinations:
            avg_consco = []
            for centroid, subtarget in zip(combination, list(all_valid_searchmodels.keys())):
                current_consco = valid_searchmodels.loc[
                    (valid_searchmodels.CENTROID_ID == centroid) & (valid_searchmodels.SUBTRGT_ID == subtarget)].CON_SCO
                avg_consco.append(float(current_consco))

            avg_consco = sum(avg_consco) / len(all_valid_searchmodels.keys())

            # Only append the combination with this score if it is higher than the one already stored
            key = ' '.join(sorted(combination))
            if key not in searchmodel_combinations.keys():
                searchmodel_combinations[key] = avg_consco
            elif avg_consco > searchmodel_combinations[key]:
                searchmodel_combinations[key] = avg_consco

        ranked_combination = sorted(searchmodel_combinations.items(), key=lambda kv: (kv[1], kv[0]))

        self.ranked_searchmodels = pd.DataFrame(ranked_combination)
        self.ranked_searchmodels.columns = ['searchmodels', 'consco']
        self.ranked_searchmodels.sort_values(by='consco', inplace=True, ascending=False)
