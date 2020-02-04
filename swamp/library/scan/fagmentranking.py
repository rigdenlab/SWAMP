import os
import swamp
import conkit.io
import itertools
import pandas as pd
from swamp.logger.swamplogger import SwampLogger
from swamp.scan.targetsplit import SplitTarget
from swamp.scan.mapalignscan import MapAlignScan
from swamp.scan.aleigenscan import AlEigenScan
from swamp.library_tools.swamplibrary import SwampLibrary
from swamp.library_tools.pdb_tools import renumber_hierarchy


class FragmentRanking(object):
    """Class to manage the rank of fragments in the library as potentially successful search models for a given set
    of predicted contacts for a target.

    Using CMO alignment tools determine the best fragments in the library to be used as search models. The target
    will be split into several subtargets (one for each helical pair with enough interhelical contact information,
    and several contact scan instances will be created.

    :param str workdir: working directory where the MR tasks will be executed
    :param str conpred: contact prediction file of the target
    :param str sspred: secondary structure prediction file of the target (must be topcons file)
    :param str conformat: format of the contact prediction file for the target
    :param str nthreads: number of parallel threads to use for CMO calculations (default: 1)
    :param list template_subset: set of templates to be used rather than the full fragment library (deafult: None)
    :param str target_pdb_benchmark: target's pdb file for benchmark purposes (default: None)
    :param str alignment_algorithm_name: algorithm used for CMO calculation (default: 'mapalign')
    :ivar bool error: True if errors have occurred at some point on the pipeline

    :example

    >>> from swamp.scan.fagmentranking import FragmentRanking
    >>> my_rank = FragmentRanking('<workdir>', '<conpred>', '<sspred>')
    >>> my_rank.rank()
    >>> my_rank.rank_searchmodels()
    """

    def __init__(self, workdir, conpred, sspred, conformat="psicov", nthreads=1, template_subset=None,
                 target_pdb_benchmark=None, alignment_algorithm_name='mapalign', logger=None):
        self._workdir = workdir
        self._conpred = conpred
        self._sspred = sspred
        self._conformat = conformat
        self._nthreads = nthreads
        self._make_workdir()
        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(logfile=os.path.join(self.workdir, "swamp_rank.log"), use_console=True,
                             console_level='info',
                             logfile_level='debug')
        else:
            self._logger = logger
        self._results = None
        self._ranked_searchmodels = None
        self._target_pdb_benchmark = target_pdb_benchmark
        self._alignment_algorithm_name = alignment_algorithm_name
        self._template_subset = template_subset
        self._splitter = None
        self._con_precision_dict = None
        self.__metadata_scanresult_raw = None

    def __repr__(self):

        return '{}(workdir="{_workdir}", conpred="{_conpred}", sspred="{_sspred}", conformat="{_conformat}", ' \
               'nthreads="{_nthreads}", template_subset="{_template_subset}", logger="{_logger}", ' \
               'target_pdb_benchmark="{_target_pdb_benchmark}", alignment_algorithm_name="{_alignment_algorithm_name}' \
               '")'.format(self.__class__.__name__, **self.__dict__)

        # ------------------ Some properties ------------------

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
    def logger(self):
        """Property to store the logger interface :obj:`~swamp.logger.swamplogger.SwampLogger`"""
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def _metadata_scanresult_raw(self):
        """Property to store the raw metadata of the scans"""
        return self.__metadata_scanresult_raw

    @_metadata_scanresult_raw.setter
    def _metadata_scanresult_raw(self, value):
        self.__metadata_scanresult_raw = value

    @property
    def con_precision_dict(self):
        """Property to store the precision of the target's predicted contacts (only used when pdb_benchmark is set)"""
        return self._con_precision_dict

    @con_precision_dict.setter
    def con_precision_dict(self, value):
        self._con_precision_dict = value

    @property
    def splitter(self):
        """Target splitter :obj:`~swamp.scan.targetsplit.SplitTarget`"""
        return self._splitter

    @splitter.setter
    def splitter(self, value):
        self._splitter = value

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
    def results(self):
        """Property to store the results obtained from the CMO scan"""
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

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
    def _tmp_cmap(self):
        return os.path.join(self.workdir, "tmp_cmap.map")

    @property
    def _tmp_pdb(self):
        if self.target_pdb_benchmark is not None:
            return os.path.join(self.workdir, 'tmp_strct.pdb')
        else:
            return None

    @property
    def _alignment_tool(self):
        """Wrapper class for the CMO alignment to be used during the ranking"""

        if self.alignment_algorithm_name == 'aleigen':
            return AlEigenScan
        elif self.alignment_algorithm_name == 'mapalign':
            return MapAlignScan
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.alignment_algorithm_name)

    @property
    def _column_reference(self):
        """Fields of the result table for each of the CMO algorithms"""

        return {
            'aleigen': ["SUBTRGT_RANK", "SUBTRGT_ID", "N_CON_MAP_A", "MAP_A", "MAP_B", "CON_SCO", "C1", "C2", "CMO",
                        "ALI_LEN", "QSCORE", "RMSD", "SEQ_ID", "N_ALIGN"],
            'mapalign': ["SUBTRGT_RANK", "SUBTRGT_ID", "N_CON_MAP_A", "MAP_A", "MAP_B", "CON_SCO", "GAP_SCO",
                         "TOTAL_SCO", "ALI_LEN", "QSCORE", "RMSD", "SEQ_ID", "N_ALIGN"]
        }

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

    # ------------------ Some general methods ------------------

    def _make_workdir(self):
        """Create the working directory"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    def _make_dataframe(self, results, **kwargs):
        """Convert the results into a :obj:`pd.Dataframe`

        :param results: Nested list with the results of each of the CMO scans
        :type results: list
        :param kwargs: Passed to :obj:`~swamp.scan.fragrank._get_best_alignment`
        :return: no value
        :rtype: None
        """

        self.results = pd.DataFrame(results)
        self.results.columns = self._column_reference[self.alignment_algorithm_name]
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

    def rank_searchmodels(self, ncontacts_threshold=28, consco_threshold=0.75, combine_searchmodels=False):
        """Get fragments groupings ranked by their combined scores. Takes in consideration same combination of fragments
        may appear more than once with fragments in different subtargets (then the combined_consco is the maximum
        possible)

        :param ncontacts_threshold: no. ofinterhelical contacts to consider the search models of subtarget (default 28)
        :type ncontacts_threshold: int
        :param consco_threshold: contact score threshold to consider a CMO alignment valid (default 0.75)
        :type consco_threshold: float
        :param combine_searchmodels: if True combine search models matching different subtargets (default False)
        :type combine_searchmodels: bool
        """

        if self.results is None:
            self.logger.error('Need to rank the fragments first!')
            return

        # Get the valid searchmodels (subtarget with enough contacts and fragment with enough consco)
        valid_searchmodels = self.results.loc[(self.results.N_CON_MAP_A >= ncontacts_threshold) & (self.results.CON_SCO >= consco_threshold),]
        if not combine_searchmodels:
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
        elif len(set(list(valid_searchmodels.SUBTRGT_ID))) == 0:
            self.logger.warning('None of the subtargets had valid searchmodels!')
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

    def rank(self, n_contacts_threshold=28):
        """Rank the fragments in the library according to their contact similarity with the different helical pairs
         in the target.

         :param n_contacts_threshold: number of interhelical contacts to consider the result CMO alignments of a subtrgt
         :type n_contacts_threshold: int
         """

        self.logger.info(self.scan_header)

        tmp_results = []
        # Get the subtargets
        self.con_precision_dict = {}
        self.logger.info("Splitting the target into sets of contacting helical pairs\n")
        self.splitter = SplitTarget(workdir=self.workdir, conpred=self.conpred, sspred=self.sspred, logger=self.logger,
                                    conformat=self.conformat, pdb_benchmark=self.target_pdb_benchmark)
        self.splitter.split()

        # Scan the library with each subtarget
        for idx, subtarget in enumerate(self.splitter.ranked_subtargets):

            if subtarget.ncontacts >= n_contacts_threshold:
                self.logger.info("Scanning library with subtarget %s (%s contacts)" % (idx + 1, subtarget.ncontacts))
                conkit.io.write(fname=self._tmp_cmap, format=self.library_format, hierarchy=subtarget)

                if self.target_pdb_benchmark is not None:
                    renumber_hierarchy(self.splitter.subtargets_pdb[subtarget.id])
                    self.splitter.subtargets_pdb[subtarget.id].write_pdb(self._tmp_pdb)
                    perfect_contacts = conkit.io.read(self._tmp_pdb, 'pdb').top_map
                    subtarget.sequence = perfect_contacts.sequence.deepcopy()
                    subtarget.set_sequence_register()
                    precision = subtarget.match(perfect_contacts).precision
                    self.con_precision_dict[subtarget.id] = precision

                scanner = self._alignment_tool(workdir=self.workdir, query=self._tmp_cmap, nthreads=self.nthreads,
                                               pdb_library=swamp.FRAG_PDB_DB, library_format=self.library_format,
                                               template_library=self.template_library, con_format=self.library_format,
                                               query_pdb_benchmark=self._tmp_pdb, template_subset=self.template_subset,
                                               logger=self.logger)
                scanner.run()

                for x in scanner.results.value:
                    tmp_results.append([idx + 1, subtarget.id, subtarget.ncontacts] + x)
                os.remove(self._tmp_cmap)

            else:
                self.logger.info("The number of interhelical contacts for subtarget %s is too small for a precise"
                                 " alignment! (%s contacts)\n" % (idx + 1, subtarget.ncontacts))

        self._metadata_scanresult_raw = tmp_results
        self._make_dataframe(tmp_results)
