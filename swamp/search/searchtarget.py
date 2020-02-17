import os
import swamp
import joblib
import conkit.io
import itertools
import pandas as pd
from pyjob import TaskFactory
from swamp.logger import SwampLogger
from swamp.search.searchjob import SearchJob
from swamp.utils.swamplibrary import SwampLibrary
from swamp.utils import TargetSplit, renumber_hierarchy


class SearchTarget(object):
    """Class to search the SWAMP library and rank search models according to their CMO with a given target.

    Using CMO alignment tools determine the best fragments in the library to be used as search models for a given \
    target. First, the target will be split into several subtargets (one for each helical pair with enough \
    interhelical contact information, and several :py:obj:`~swamp.search.searchjob.SearchJob` instances will be \
    executed.

    :param str workdir: working directory for the :py:obj:`~swamp.search.searchjob.SearchJob` instances
    :param str conpred: contact prediction file of the target
    :param str sspred: secondary structure prediction file of the target (must be topcons format file)
    :param str conformat: format of the contact prediction file provided for the target (default: 'psicov')
    :param str nthreads: number of parallel threads to use in the library search (default: 1)
    :param tuple template_subset: set of templates to be used instead of the full fragment library (deafult: None)
    :param str target_pdb_benchmark: provide a target's pdb file for benchmark purposes (default: None)
    :param str alignment_algorithm_name: algorithm used for CMO calculation (default: 'mapalign')
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the search (default None)
    :param int n_contacts_threshold: min. no. of interhelical cont. to include a subtarget in the search (default: 28)
    :param str platform: scheduler system where the array will be executed (default 'sge')
    :param str queue_name: name of the scheduler queue where the tasks should be sent (default None)
    :param str queue_environment: name of the scheduler environment where the tasks should be sent (default None)
    :param str python_interpreter: python interpreter to be used for task execution (default '$CCP4/bin/ccp4-python')
    :ivar str shell_interpreter: shell interpreter to be used for task execution (default '/bin/bash')
    :ivar bool error: True if errors have occurred at some point in the pipeline
    :ivar :py:obj:`~swamp.utils.targetsplit.TargetSplit` target: contains information about the target and subtargets
    :ivar dict con_precision_dict: a dictionary with the contact precission for each given subtarget prediction
    :ivar dict search_pickle_dict: a dictionary with the :py:attr:`~swamp.search.searchjob.SearchJob.pickle_fname` \
     created by each :py:obj:`~swamp.search.searchjob.SearchJob` instance in this search
    :ivar list results: a nested list with the results obtained in the search againts the library
    :ivar list scripts: a list with the instances of :py:obj:`pyjob.Script` that will be executed to complete the search
    :ivar :py:obj:`pandas.DataFrame` ranked_searchmodels: a dataframe with the search models ranked by the CMO \
    obtained in the search

    :example:

    >>> from swamp.search import SearchTarget
    >>> my_rank = SearchTarget('<workdir>', '<conpred>', '<sspred>')
    >>> my_rank.search()
    >>> my_rank.rank()
    """

    def __init__(self, workdir, conpred, sspred, conformat="psicov", nthreads=1, template_subset=None, logger=None,
                 target_pdb_benchmark=None, alignment_algorithm_name='mapalign', n_contacts_threshold=28,
                 python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'),
                 platform='sge', queue_name=None, queue_environment=None):
        self.workdir = workdir
        self.conpred = conpred
        self.sspred = sspred
        self.conformat = conformat
        self.nthreads = nthreads
        self.target_pdb_benchmark = target_pdb_benchmark
        self.alignment_algorithm_name = alignment_algorithm_name
        self.template_subset = template_subset
        self.n_contacts_threshold = n_contacts_threshold
        self.platform = platform
        self.queue_name = queue_name
        self.queue_environment = queue_environment
        self.python_interpreter = python_interpreter
        self.shell_interpreter = '/bin/bash'
        self.con_precision_dict = None
        self.search_pickle_dict = None
        self.scripts = None
        self.results = None
        self.ranked_searchmodels = None
        self._make_workdir()

        if logger is None:
            self.logger = SwampLogger(__name__)
            self.logger.init(logfile=os.path.join(self.workdir, "swamp_rank.log"), use_console=True,
                             console_level='info',
                             logfile_level='debug')
        else:
            self.logger = logger
        self.target = TargetSplit(workdir=self.workdir, conpred=self.conpred, sspred=self.sspred, logger=self.logger,
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
    def search_header(self):
        """Header displayed when initiating :py:obj:`~swamp.logger.swamplogger.SwampLogger`"""

        return """**********************************************************************
*****************            SWAMP SEARCH            *****************
**********************************************************************

"""

    @property
    def _tmp_cmap(self):
        """A temporary file name to contain the contact predicionts of each subtarget at \
        :py:attr:`swamp.utils.targetsplit.TargetSplit.ranked_subtargets`"""
        return os.path.join(self.workdir, "tmp_cmap_{}.map")

    @property
    def _search_workdir(self):
        """The workind girectory for each :py:obj:`~swamp.search.searchjob.SearchJob` instance in this search"""
        return os.path.join(self.workdir, "search_{}")

    @property
    def _tmp_pdb(self):
        """A temporary file name to contain the pdb file each subtarget at \
        :py:attr:`swamp.utils.targetsplit.TargetSplit.ranked_subtargets`"""
        if self.target_pdb_benchmark is not None:
            return os.path.join(self.workdir, 'tmp_strct_{}.pdb')
        else:
            return None

    @property
    def template_library(self):
        """Location of the template library to be used with the :py:obj:`~swamp.search.searchjob.SearchJob` instances"""

        if self.alignment_algorithm_name == 'aleigen':
            return swamp.FRAG_ALEIGEN_DB
        elif self.alignment_algorithm_name == 'mapalign':
            return swamp.FRAG_MAPALIGN_DB
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    @property
    def library_format(self):
        """Dictionary specifying the template library format to be used in each \
        :py:obj:`~swamp.search.searchjob.SearchJob` instance"""

        if self.alignment_algorithm_name == 'aleigen':
            return 'aleigen'
        elif self.alignment_algorithm_name == 'mapalign':
            return 'mapalign'
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    @property
    def _other_task_info(self):
        """Dictionary with the extra **kwags passed to :py:obj:`pyjob.TaskFactory`"""

        info = {'directory': self.workdir, 'shell': self.shell_interpreter, 'name': 'swamp_search'}

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
        """A list of column names for :py:attr:`~swamp.search.searchtarget.SearchTarget.results`"""

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
        """Create the :py:attr:`~swamp.search.searchtarget.SearchTarget.workdir`"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    def _create_scripts(self):
        """Populate the :py:attr:`~swamp.search.searchtarget.SearchTarget.scripts` list with all the \
        :py:obj:`pyjob.Script` instances that need to be executed. There is one instance for each of the subtargets \
        at :py:attr:`swamp.utils.targetsplit.TargetSplit.ranked_subtargets` that passes the selected \
        :py:attr:`~swamp.search.searchtarget.SearchTarget.n_contacts_threshold` threshold"""

        self.scripts = []
        self.con_precision_dict = {}
        self.search_pickle_dict = {}

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

                searcher = SearchJob(**self._search_info(idx))

                self.scripts.append(searcher.script)
                self.search_pickle_dict[searcher.pickle_fname] = subtarget

            else:
                self.logger.info("No. of contacts below the %s threshold. "
                                 "Subtarget will not be used in the search." % self.n_contacts_threshold)

        self.logger.info('%s subtargets will be used in the search. Creating task now.' % len(self.search_pickle_dict))

    def _search_info(self, idx):
        """Create **kwargs passed to a :py:obj:`~swamp.search.searchjob.SearchJob` instance

        :param int idx: the index of the search job, used as :py:attr:`~swamp.search.searchjob.SearchJob.id`
        :returns: dictionary with the **kwargs to be passed into :py:obj:`~swamp.search.searchjob.SearchJob` instance
        """

        info = {'workdir': self._search_workdir.format(idx), 'pdb_library': swamp.FRAG_PDB_DB,
                'query': self._tmp_cmap.format(idx), 'algorithm': self.alignment_algorithm_name,
                'template_subset': self.template_subset, 'python_interpreter': self.python_interpreter,
                'template_library': self.template_library, 'library_format': self.library_format,
                'logger': self.logger, 'con_format': self.library_format, 'id': idx}

        if self.target_pdb_benchmark is not None:
            info['query_pdb_benchmark'] = self._tmp_pdb.format(idx),

        return info

    def _make_dataframe(self, results, **kwargs):
        """Convert the :py:attr:`~swamp.search.searchtarget.SearchTarget.results` into a :py:obj:`pandas.Dataframe`

        :param list results: Nested list with the results of each of the CMO scans
        :param kwargs: arguments are passed to :py:func:`~swamp.search.seachtarget.SearchTarget._get_best_alignment`
        """

        self.results = pd.DataFrame(results)
        self.results.columns = self._column_reference
        self.results["CENTROID_ID"] = [os.path.basename(fname).split('.')[0] for fname in self.results.MAP_B.to_list()]
        self._get_best_alignment(**kwargs)

    def _get_best_alignment(self, select_by=("CENTROID_ID", "SUBTRGT_ID")):
        """Update the :py:attr:`~swamp.search.searchtarget.SearchTarget.results` dataframe to include only the results \
        with the optimal CMO alignment between helical pairs. This method considers inverted fragments and takes \
        highest CMO score to determine the optimal alignment

         :param list select_by: indicate the figure of merit by which the alignments will be grouped
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
        """Recover the results from all the :py:attr:`~swamp.search.searchjob.SearchJob.pickle_fname` indicated in \
        :py:attr:`~swamp.search.searchtarget.SearchTargert.search_pickle_dict`

        :returns a list with the results loaded from the pickle files at \
        :py:attr:`~swamp.search.searchtarget.SearchTargert.search_pickle_dict`
        """

        results = []
        for pickle_fname, subtarget in zip(self.search_pickle_dict.keys(), self.search_pickle_dict.values()):

            if os.path.isfile(pickle_fname):
                self.logger.debug('Retrieving results from %s' % pickle_fname)
                current_job_results = joblib.load(pickle_fname)
                for result in current_job_results:
                    results.append([self.target.ranked_subtargets.index(subtarget) + 1, subtarget.id,
                                    subtarget.ncontacts] + result)

            else:
                self.logger.warning('Cannot find pickle file! %s' % pickle_fname)

        return results

    def search(self):
        """Search the library by calculating the CMO between the observed contacts and the  query predicted contacts

        This method will run a :py:obj:`~swamp.search.searchjob.SearchJob` instance for each subtarget \
        at :py:attr:`swamp.utils.targetsplit.TargetSplit.ranked_subtargets` that passes the selected \
        :py:attr:`~swamp.search.searchtarget.SearchTarget.n_contacts_threshold` threshold
        """

        self.logger.info(self.search_header)
        self.logger.info("Splitting the target into sets of contacting helical pairs")
        self.target.split()
        if self.target.error:
            self.logger.warning('Previous errors prevent scanning the library with the target contacts!')
            return

        self.logger.info('Creating a list of jobs to search the library using contacts.')
        self._create_scripts()
        self.logger.info('Sending jobs now.')

        with TaskFactory(self.platform, tuple(self.scripts), **self._other_task_info) as task:
            task.run()
            self.logger.info('Waiting for workers...')

        self.logger.info('All search tasks have been completed! Retrieving results')
        results = self.recover_results()
        self._make_dataframe(results)

    def rank(self, consco_threshold=0.75, combine_searchmodels=False):
        """ Get the top search models as indicated by the results of the search. This method can also try to combine \
        multiple search models by adding their CMOs. Takes in consideration same combination of fragments may appear \
        more than once with fragments in matching subtargets.

         :param float consco_threshold: CMO threshold to consider an alignment valid (default 0.75)
         :param bool combine_searchmodels: if True combine search models matching different subtargets (default False)
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
