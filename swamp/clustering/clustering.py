import os
import abc
import swamp
import threading
import numpy as np
import pandas as pd
from statistics import mean
import swamp.utils.swamplibrary
from swamp.logger import SwampLogger
from swamp.wrappers.gesamt import Gesamt
from itertools import groupby, combinations
from sklearn.metrics import silhouette_score
from sklearn.model_selection import ParameterSampler

ABC = abc.ABCMeta('ABC', (object,), {})


class Clustering(ABC):
    """Abstract class for SWAMP library clustering purposes, used to cluster similar helical pairs into ensembles

    Implements data structures and methods commonly used in clustering tasks and optimal hyper-parameter search

    :param int nthreads: number of threads to be used for fragment clustering (default 1)
    :param int n_iter: number of iterations for the :py:func:`~swamp.clustering.clustering.grid_search` of the optimal \
     clustering hyper paramaters (default 100)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logger instance to record log messages
    :ivar bool error: True if errors have occurred at some point on the process of clustering
    :ivar dict best_params: contains the optimal hyper-parameters as found on the \
    :py:func:`~swamp.clustering.clustering.grid_search`
    :ivar list labels: a list with the cluster labels assigned to each member of the \
    :py:attr:`~swamp.clustering.clustering.similarity_mtx`
    :ivar list child_threads: a list with the :py:attr:`threading.thread.name` of each child \
    :py:obj:`threading.thread` instance used on the :py:func:`~swamp.clustering.clustering.grid_search`
    :ivar `pandas.DataFrame` rmsd_mtx: square dataframe with the rmsd distance across framgents in the library
    :ivar `pandas.DataFrame` similarity_mtx: square dataframe with the similarity across framgents in the library
    :ivar `pandas.DataFrame` nalign_mtx: square dataframe with the no. of aligned residues between framgents in the \
    library
    :ivar dict centroid_dict: dictionary with the centroids associated with each cluster id
    :ivar dict cluster_dict: dictionary with the list of fragments ids that form each cluster
    :ivar dict composition_dict: dictionary with the final composition of each ensemble
    :ivar dict ensemble_dict: dictionary with the list fragments ids that form each ensemble
    :ivar `threading.Semaphore` semaphore: a `threading.Semaphore` instance to control the execution of threads \
    on :py:func:`~swamp.clustering.clustering.grid_search`
    :ivar `~swamp.clustering.clustering.ParameterSearchResults` gridsearch_results: a \
    :py:obj:`~swamp.clustering.clustering.ParameterSearchResults` instance with the results obtained on \
    :py:func:`~swamp.clustering.clustering.grid_search`
    """

    def __init__(self, nthreads=1, n_iter=100, logger=None):
        self.error = False
        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(use_console=True, console_level='info')
        else:
            self._logger = logger
        self.best_params = None
        self.labels = None
        self.child_threads = None
        self.nthreads = nthreads
        self.n_iter = n_iter
        self.cluster_dict = None
        self.centroid_dict = None
        self.ensemble_dict = None
        self.rmsd_mtx = None
        self.nalign_mtx = None
        self.similarity_mtx = None
        self.semaphore = threading.Semaphore(value=self.nthreads)
        self.gridsearch_results = ParameterSearchResults(self.logger)
        self.composition_dict = {}

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def cluster(self):
        """ Abstract method to run the clustering algorithm"""
        pass

    @abc.abstractmethod
    def _clustering(self, **kwargs):
        """Abstract method to hold the clustering algorithm"""
        pass

    @property
    @abc.abstractmethod
    def _algorithm_name(self):
        """Abstract property to hold the name of the specific clustering algorithm"""
        pass

    @property
    @abc.abstractmethod
    def _hyper_params(self):
        """Abstract property to hold the hyperparameters for the specific clustering algorithm"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def clustering_header(self):
        """The header to be displayed when initialising the logger"""

        return """**********************************************************************
*****************         SWAMP CLUSTERING           *****************
**********************************************************************

"""

    @property
    def _param_sample(self):
        """Sample of hyperparameters to be tested in the :py:func:`~swamp.clustering.clustering.grid_search`"""

        sampler = ParameterSampler(self._hyper_params, n_iter=self.n_iter, random_state=41)
        return [d for d in sampler]

    @property
    def distance_mtx(self):
        """Square matrix that corresponds to 1 - :py:attr:`~swamp.clustering.clustering.distance_mtx`"""
        return 1 - self.similarity_mtx

    @property
    def fragment_list(self):
        """Columns of the :py:attr:`~swamp.clustering.clustering.distance_mtx`, which correspond with the list of \
        fragments in the library"""
        return self.distance_mtx.columns

    @property
    def _reference_input(self):
        """Dictionary with the ihe input required for each clustering algorithm"""

        return {"affinity": self.similarity_mtx,
                "spectral": self.similarity_mtx,
                "optics": self.distance_mtx,
                "agglomerative": self.distance_mtx,
                "dbscan": self.distance_mtx
                }

    # ------------------ Some protected methods ------------------

    def _join_threads(self):
        """Join the threads listed at :py:attr:`~swamp.clustering.child_threads` to wait for their completion"""

        mainthread = threading.current_thread()
        self.logger.debug("Mainthread is %s" % mainthread.getName())
        for t in threading.enumerate():
            if t is not mainthread and t.getName() in self.child_threads:
                self.logger.debug("joining thread %s" % t.getName())
                t.join()

    def _get_clst_qscore(self, frag_ids):
        """Method to get the average qscore among the models in a given cluster

        :param tuple frag_ids: a list with the fragment ids of contained in the cluster
        :returns: the average qscore among the models in the cluster (float)
        """

        qscores = []
        for frag_pair in combinations(frag_ids, 2):
            frag_a_scores = self.similarity_mtx[frag_pair[0]]
            qscores.append(frag_a_scores[frag_a_scores.index == frag_pair[1]].values[0])
        return mean(qscores)

    def _test_params(self, idx, **kwargs):
        """Method to test clustering the library of fragments with a given set of parameters

        :param int idx: index of the current test within the :py:func:`~swamp.clustering.clustering.grid_search`
        :param dict kwargs: arguments passed to :py:obj:`~swamp.cluster.clustering._clustering()`
        """

        self.semaphore.acquire()
        self.logger.debug("Grid search parameters %s" % dict(**kwargs))
        clst = self._clustering(**kwargs)
        clst.fit(self._reference_input[self._algorithm_name])
        self.logger.debug("Model fit done %s" % dict(**kwargs))
        self._gridsearch_results.register(
            new_results=pd.DataFrame([[idx] + self.assess_clustering(labels=clst.labels_)],
                                     columns=["params_idx", "n_clst", "clst_size", "clst_qscore", "silhouette_score",
                                              "singletons"]))
        self.logger.debug("Done with grid search %s" % dict(**kwargs))
        self.semaphore.release()

    # ------------------ Some general methods ------------------

    def get_clst_id(self, frag_id):
        """Get the unique cluster identifier where a given fragment is assigned

        :param str frag_id: the id of the fragment of interest
        :returns: the cluster id where the fragment can be found (str)
        """

        if self.cluster_dict is None:
            return None

        for clst_id in self.cluster_dict.keys():
            if frag_id in self.cluster_dict[clst_id]:
                return clst_id

    def get_centroid_id(self, frag_id):
        """Get the unique cluster identifier where a given centroid is assigned

        :param str frag_id: the id of the centroid of interest
        :returns: the cluster id where the fragment can be found (str)
        """

        if self.centroid_dict is None:
            return None

        for clst_id in self.centroid_dict.keys():
            if frag_id == self.centroid_dict[clst_id]:
                return clst_id

    def register_library(self, library):
        """Register a given :py:obj:`~swamp.utils.swamplibrary.SwampLibrary` instance in order to load the fragment \
        distance info

        :argument `~swamp.utils.swamplibrary.SwampLibrary` library: the \
        :py:obj:`~swamp.utils.swamplibrary.SwampLibrary` insntance to be registered
        :raises TypeError: if `library` is not a :py:obj:`~swamp.utils.swamplibrary.SwampLibrary` insntance
        """

        if not isinstance(library, swamp.utils.swamplibrary.SwampLibrary):
            raise TypeError('Library to be registered must be a swamp.utils.swamplibrary.SwampLibrary !')

        nalign_mtx = library.nalign_matrix.fillna(0, inplace=False)
        self.nalign_mtx = nalign_mtx.astype(int)

        rmsd_mtx = library.rmsd_matrix.fillna(0, inplace=False)
        self.rmsd_mtx = rmsd_mtx.astype(float)

        qscore_mtx = library.qscore_matrix.fillna(0, inplace=False)
        self.similarity_mtx = qscore_mtx.astype(float)

    def grid_search(self):
        """Method to do a grid random search for a range of clustering hyper-parameters defined at \
        :py:attr:`~swamp.clustering.clustering._hyper_params`"""

        self.logger.info("Starting grid hyper parameter search with %s parallel threads and %s iterations" % (
            self.nthreads, self.n_iter))
        self.child_threads = []
        for idx, params in enumerate(self._param_sample):
            t = threading.Thread(target=self._test_params, args=(idx,), kwargs=(params))
            self.logger.debug("Sending thread %s" % t.getName())
            self.child_threads.append(t.getName())
            t.start()
        self.logger.info("Waiting for workers")
        self._join_threads()
        self.logger.debug("Grid parameter search is done!")

    def assess_clustering(self, labels):
        """Method to asses the quality of the obtained clustering

        :param tuple labels: the labels assigned to each fragment
        :returns: a tuple with no. clusters, average cluster size, average cluster qscore, silhouette score and the \
         no. singletons
        """

        # Clst size
        tmp_labels = [x for x in labels if x != -1]
        if len(tmp_labels) == 0:
            return [np.nan] * 5
        tmp_labels.sort()
        clst_size = mean([len(list(group)) for key, group in groupby(tmp_labels)])
        # Singletons
        singletons = list(labels).count(-1)
        # Clst qscore
        clst_dict = self.get_cluster_dict(labels=labels, inplace=False)
        qscore_list = []
        for clst in clst_dict.keys():
            if clst != -1:
                qscore_list.append(self._get_clst_qscore(clst_dict[clst]))
        clst_qscore = mean(qscore_list)
        # Silhouette score
        silhouette_sco = silhouette_score(self.distance_mtx, metric='precomputed', labels=labels)

        return [len(set([x for x in tmp_labels])), clst_size, clst_qscore, silhouette_sco, singletons]

    def get_cluster_dict(self, labels, inplace=True):
        """Method to generate a cluster dictionary containing the cluster identifiers as keys and the fragment names \
        as the values

        :param tuple labels: the labels assigned to each fragment in the library
        :param bool inplace: if True, then it sets :py:attr:`~swamp.clustering.clustering.cluster_dict` (default True)
        :returns: a dictionary containing the cluster id as keys and the frag ids as the values (if not `inplace`)
        """

        self.logger.debug('Creating cluster dictionary')

        clst_dict = {}
        for idx, label in enumerate(labels):
            if label in clst_dict.keys():
                clst_dict[label].append(self.similarity_mtx.columns[idx])
            else:
                clst_dict[label] = [self.similarity_mtx.columns[idx]]

        if inplace:
            self.cluster_dict = clst_dict
        else:
            return clst_dict

    def get_centroid_dict(self):
        """Get a dictionary with the centroids of each cluster at \
        :py:attr:`~swamp.clustering.clustering.cluster dict`"""

        self.logger.debug('Creating centroid dictionary')

        if self.cluster_dict is None:
            return None

        self.centroid_dict = {}
        for clst_id in self.cluster_dict.keys():
            if clst_id != -1 and len(self.cluster_dict[clst_id]) > 1:
                if len(self.cluster_dict[clst_id]) == 2:
                    self.centroid_dict[clst_id] = self.cluster_dict[clst_id][0]
                else:
                    current_centroid = (0, None)
                    # Get the average qscore respect the other fragments in the cluster (except the current candidate)
                    for fragment in self.cluster_dict[clst_id]:
                        tmp_mtx = self.similarity_mtx[fragment].drop(fragment, 0)
                        tmp_mtx = tmp_mtx.drop([x for x in tmp_mtx.index.values if x not in self.cluster_dict[clst_id]],
                                               0)
                        avg_qscore = mean(list(tmp_mtx))
                        # The candidate with the lowest avg qscore is the centroid
                        if avg_qscore > current_centroid[0]:
                            current_centroid = (avg_qscore, fragment)
                    self.centroid_dict[clst_id] = current_centroid[1]

    def get_ensemble_dict(self, rmsd_threshold=0.7, nalign_threshold=30, qscore_threshold=None, nthreads=1):
        """Merge similar fragments to create a dictionary with fragment ensembles (clustering with replacement)

        An rmsd and no. aligned residues threshold can be used, but if a qscore threshold is given, this one will
        be used instead

        :param float rmsd_threshold: the rmsd threshold at which fragments are included into the ensemble (default 0.7)
        :param int nalign_threshold: threshold of aligned residues to include a fragment into an ensemble (default 30)
        :param float qscore_threshold: qscore threshold at which fragments are included into the ensemble (default None)
        :param int nthreads: number of threads to compute the ensemble optimal parameter (default 1)
        """

        self.logger.debug('Creating ensemble dictionary')
        self.ensemble_dict = {}

        # Check all centroids
        for clst_id in self.centroid_dict.keys():
            if clst_id != -1:
                self.ensemble_dict[clst_id] = self.get_ensemble(centroid_id=clst_id, rmsd_threshold=rmsd_threshold,
                                                                nalign_threshold=nalign_threshold, nthreads=nthreads,
                                                                qscore_threshold=qscore_threshold)

    def get_ensemble(self, centroid_id, rmsd_threshold=0.7, nalign_threshold=30, qscore_threshold=None, nthreads=1):
        """Search for fragments to form an ensemble given a centroid identifier

        Fragments within the threshold distance from the centroid will be included in the ensemble. Other centroids \
        will be excluded from this search.

        :param int centroid_id: centroid identifier of the centroid of interest
        :param float rmsd_threshold: the rmsd threshold at which fragments are included into the ensemble (default 0.7)
        :param int nalign_threshold: threshold of aligned residues to include a fragment into an ensemble (default 30)
        :param float qscore_threshold: qscore threshold at which fragments are included into the ensemble (default None)
        :param int nthreads: number of threads to compute the ensemble optimal model alignment (default 1)
        :returns: a :py:obj:`gemmi.Structure` hierarchy with the ensemble
        """

        if self.centroid_dict is None:
            self.get_centroid_dict()

        def _lower_threshold(threshold, logger):
            threshold -= 0.03
            logger.debug('No fragments found, lowering qscore threshold to %s' % threshold)
            return list(set(tmp_dist_mtx.index[tmp_dist_mtx.astype(float) >= threshold].tolist())), threshold

        def _increase_threshold(threshold, logger):
            threshold += 0.03
            logger.debug('Found %s fragments, increasing qscore threshold to %s' % (
                len(cluster_fragments), threshold))
            return list(set(tmp_dist_mtx.index[tmp_dist_mtx.astype(float) >= threshold].tolist())), threshold

        self.logger.debug('Retrieving ensemble for cluster %s' % centroid_id)
        fname = os.path.join(swamp.FRAG_PDB_DB, '{}.pdb')
        centroid = self.centroid_dict[centroid_id]

        # Get the the valid fragments close enough to this centroid
        if qscore_threshold is None:
            tmp_dist_mtx = self.rmsd_mtx[centroid].drop(centroid, 0)
            good_distance_fragments = tmp_dist_mtx.index[tmp_dist_mtx.astype(float) <= rmsd_threshold].tolist()
            tmp_nalign_mtx = self.nalign_mtx[centroid].drop(centroid, 0)
            good_nalign_fragments = tmp_nalign_mtx.index[tmp_nalign_mtx.astype(int) >= nalign_threshold].tolist()
            cluster_fragments = list(set(good_nalign_fragments) & set(good_distance_fragments))
        else:
            tmp_dist_mtx = self.similarity_mtx[centroid].drop(centroid, 0)
            good_qscore_fragments = tmp_dist_mtx.index[tmp_dist_mtx.astype(float) >= qscore_threshold].tolist()
            cluster_fragments = list(set(good_qscore_fragments))
            # This optimizes results...
            while len(cluster_fragments) == 0 and qscore_threshold > 0.8:
                cluster_fragments, qscore_threshold = _lower_threshold(qscore_threshold, self.logger)
            while len(cluster_fragments) > 6 and qscore_threshold < 0.92:
                cluster_fragments, qscore_threshold = _increase_threshold(qscore_threshold, self.logger)

        # Add the centroid to the ensemble and exclude any other centroid
        cluster_fragments = list(set(cluster_fragments) - set(self.centroid_dict.values()))
        cluster_fragments.append(centroid)
        # Only add the ensemble if this combination of fragments didn't appear yet
        if sorted(cluster_fragments) not in self.composition_dict.values():
            self.composition_dict[centroid_id] = sorted(cluster_fragments)

        # If there are enough good fragments, make an ensemble and write into a file
        self.logger.debug('Retrieving optimum alignment for %s fragments [%s]' % (
            len(cluster_fragments), ', '.join(cluster_fragments)))
        gesamt, hierarchy = Gesamt.get_optimum_alignment([fname.format(frag_id) for frag_id in cluster_fragments],
                                                         nthreads=nthreads, logger=self.logger)
        if gesamt:
            self.logger.debug('Optimal qscore: %s' % gesamt.qscore)
        return hierarchy


class ParameterSearchResults(object):
    """Class to hold the results from a multi-threaded :py:func:`~swamp.clustering.clustering.grid_search`

    Implements :py:obj:`~threading.Semaphore` methods to regulate thread I/O into a result list

    :param `~swamp.logger.swamplogger.SwampLogger` logger: logger instance to record thread log messages
    :ivar `threading.lock` lock: lock to control I/O to the result instance
    :ivar `pandas.DataFrame` value: dataframe with the results of the grid search
    """

    def __init__(self, logger):
        self.lock = threading.Lock()
        self.value = pd.DataFrame(
            columns=["params_idx", "n_clst", "clst_size", "clst_qscore", "silhouette_score", "singletons"])
        self.logger = logger

    def register(self, new_results):
        """Register a given set of new results into :py:attr:`~swamp.clustering.ParameterSearchResults.value`

        :param `pandas.DataFrame` new_results: the set of new results to be registered
        """

        self.logger.debug('Waiting for lock')
        self.lock.acquire()
        self.logger.debug('Acquired lock')
        self.value = pd.concat([self.value, new_results], 0)
        self.lock.release()
        self.logger.debug('Released lock')
