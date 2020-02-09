import os
import abc
import swamp
import threading
import numpy as np
import pandas as pd
from statistics import mean
from swamp.wrappers.gesamt import Gesamt
from itertools import groupby, combinations
from sklearn.metrics import silhouette_score
from swamp.logger.swamplogger import SwampLogger
from sklearn.model_selection import ParameterSampler

ABC = abc.ABCMeta('ABC', (object,), {})


class Clustering(ABC):
    """Abstract class for SWAMP library clustering purposes

    Implements data structures and methods commonly used in clustering tasks and optimal hyperparameter search

    :param int nthreads: number of threads to be used for fragment clustering (default 1)
    :param int n_iter: number of iterations for the random grid search of the optimal clustering hyper paramaters (default 100)
    :param :obj:`swamp.logger.swamplogger` logger: logger instance to record log messages
    :ivar bool error: True if errors have occurred at some point on the pipeline
    :ivar dict best_params: contains the optimal hyperparameters for clustering
    :ivar :obj:`pandas.DataFrame` rmsd_mtx: square dataframe with the rmsd distance across framgents in the library
    :ivar :obj:`pandas.DataFrame` similarity_mtx: square dataframe with the similarity across framgents in the library
    :ivar :obj:`pandas.DataFrame` nalign_mtx: square dataframe with the no. of aligned residues between framgents in the library
    :ivar dict centroid_dict: dictionary with the centroids associated with each cluster id
    :ivar dict cluster_dict: dictionary with the fragments that form each cluster
    :ivar dict ensemble_dict: dictionary with the fragments that form each ensemble
    :ivar bool error: True if errors have occurred at some point on the pipeline
    """

    def __init__(self, nthreads=1, n_iter=100, logger=None):
        self._error = False
        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(use_console=True, console_level='info')
        else:
            self._logger = logger
        self._best_params = None
        self._labels = None
        self._child_threads = None
        self._nthreads = nthreads
        self._n_iter = n_iter
        self._outdir = None
        self._cluster_dict = None
        self._centroid_dict = None
        self._ensemble_dict = None
        self._rmsd_mtx = None
        self._nalign_mtx = None
        self._similarity_mtx = None
        self._semaphore = threading.Semaphore(value=self.nthreads)
        self._gridsearch_results = ParameterSearchResults(self.logger)
        self._composition_dict = {}

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def cluster(self):
        """ Abstract method to run the wrapper"""
        pass

    @abc.abstractmethod
    def _clustering(self, **kwargs):
        """Abstract method to hold the clustering algorithm"""
        pass

    @property
    @abc.abstractmethod
    def _algorithm_name(self):
        """Abstract method to hold the name of the specific clustering algorithm"""
        pass

    @property
    @abc.abstractmethod
    def _hyper_params(self):
        """Abstract method to hold the name of the hyperparameters for the specific clustering algorithm"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def clustering_header(self):
        """Abstract property to store the clustering header for the logger"""

        return """**********************************************************************
***************         SWAMP CLUSTERING %s          *************
**********************************************************************

""" % self._algorithm_name.upper()

    @property
    def composition_dict(self):
        """Property to store the composition of the final ensembles (to output to the actual library)"""
        return self._composition_dict

    @composition_dict.setter
    def composition_dict(self, value):
        self._composition_dict = value

    @property
    def similarity_mtx(self):
        return self._similarity_mtx

    @similarity_mtx.setter
    def similarity_mtx(self, value):
        self._similarity_mtx = value

    @property
    def cluster_dict(self):
        return self._cluster_dict

    @cluster_dict.setter
    def cluster_dict(self, value):
        self._cluster_dict = value

    @property
    def ensemble_dict(self):
        return self._ensemble_dict

    @ensemble_dict.setter
    def ensemble_dict(self, value):
        self._ensemble_dict = value

    @property
    def centroid_dict(self):
        return self._centroid_dict

    @centroid_dict.setter
    def centroid_dict(self, value):
        self._centroid_dict = value

    @property
    def child_threads(self):
        return self._child_threads

    @child_threads.setter
    def child_threads(self, value):
        self._child_threads = value

    @property
    def nthreads(self):
        return self._nthreads

    @nthreads.setter
    def nthreads(self, value):
        self._nthreads = value

    @property
    def n_iter(self):
        return self._n_iter

    @n_iter.setter
    def n_iter(self, value):
        self._n_iter = value

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, value):
        self._outdir = value

    @property
    def best_params(self):
        return self._best_params

    @best_params.setter
    def best_params(self, value):
        self._best_params = value

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def semaphore(self):
        return self._semaphore

    @semaphore.setter
    def semaphore(self, value):
        self._semaphore = value

    @property
    def labels(self):
        return self._labels

    @labels.setter
    def labels(self, value):
        self._labels = value

    @property
    def nalign_mtx(self):
        return self._nalign_mtx

    @nalign_mtx.setter
    def nalign_mtx(self, value):
        self._nalign_mtx = value

    @property
    def rmsd_mtx(self):
        return self._rmsd_mtx

    @rmsd_mtx.setter
    def rmsd_mtx(self, value):
        self._rmsd_mtx = value

    @property
    def _param_sample(self):
        """Sample of hyperparameters to be tested in the randomised grid search"""

        sampler = ParameterSampler(self._hyper_params, n_iter=self.n_iter, random_state=41)
        return [d for d in sampler]

    @property
    def distance_mtx(self):
        """1 - the distances between fragments"""
        return 1 - self.similarity_mtx

    @property
    def fragment_list(self):
        """All the fragments loaded in the distance/similarity matrices"""
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
        """Join the child theads so that the code waits for their completion"""

        mainthread = threading.current_thread()
        self.logger.debug("Mainthread is %s" % mainthread.getName())
        for t in threading.enumerate():
            if t is not mainthread and t.getName() in self.child_threads:
                self.logger.debug("joining thread %s" % t.getName())
                t.join()

    def _get_clst_qscore(self, frag_ids):
        """Method to get the average qscore among the models in a given cluster

        :param frag_ids: a list with the fragment ids of contained in the cluster
        :type frag_ids: list, tuple
        :returns qscore: the average qscore among the models in the cluster
        :rtype float
        """

        qscores = []
        for frag_pair in combinations(frag_ids, 2):
            frag_a_scores = self.similarity_mtx[frag_pair[0]]
            qscores.append(frag_a_scores[frag_a_scores.index == frag_pair[1]].values[0])
        return mean(qscores)

    def _test_params(self, idx, **kwargs):
        """Method to cluster the fragments with a given set of parameters

        :param idx: index of the current test respect the rest of runs in the randomised grid search
        :type idx: int
        :param kwargs: arguments passed to the clustering instance :obj:`swamp.cluster.clustering._clustering()`
        :type kwargs: dict
        :returns nothing
        :rtype None
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
        """Get the cluster id where a fragment is located

        :param frag_id: the id of the fragment of interest
        :type frag_id: str
        :returns clst_id: the cluster id where the fragment can be found
        :rtype str, int, None
        """

        if self.cluster_dict is None:
            return None

        for clst_id in self.cluster_dict.keys():
            if frag_id in self.cluster_dict[clst_id]:
                return clst_id

    def get_centroid_id(self, frag_id):
        """Get the cluster id where of a given centroid is located

        :param frag_id: the id of the centroid of interest
        :type frag_id: str
        :returns clst_id: the cluster id where the fragment can be found
        :rtype str, int, None
        """

        if self.centroid_dict is None:
            return None

        for clst_id in self.centroid_dict.keys():
            if frag_id == self.centroid_dict[clst_id]:
                return clst_id

    def make_outdir(self):
        """Method to crete the workdir for the wrapper"""

        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)

    def register_library(self, library):
        """Register a given fragment library to load the distance matrices

        :argument library: fragment library with the distance matrices loaded
        :type :object ~swamp.framgent_library.swamplibrary.FragmentLibrary
        """

        nalign_mtx = library.nalign_matrix.fillna(0, inplace=False)
        self.nalign_mtx = nalign_mtx.astype(int)

        rmsd_mtx = library.rmsd_matrix.fillna(0, inplace=False)
        self.rmsd_mtx = rmsd_mtx.astype(float)

        qscore_mtx = library.qscore_matrix.fillna(0, inplace=False)
        self.similarity_mtx = qscore_mtx.astype(float)

    def parameter_search(self):
        """Method to do a grid random search for optimal clustering hyperparameters"""

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
        """Method to asses the quality of the clustering

        :param labels: the list containing the labels assigned to each fragment
        :type labels: list, tuple
        :returns no. clusters, average cluster size, average cluster qscore, silhouette score, no. singletons
        :rtype tuple
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
        """Method to generate a cluster dictionary containing the cluster id as keys and the frag ids as the values

        :param labels: list with the labels assigned to each fragment in the library
        :type labels: list, tuple
        :param inplace: if True, then it set the cluster dicionary attribute of the instance with the result
        :type inplace: bool
        :returns if inplace is True, a dictionary containing the cluster id as keys and the frag ids as the values
        :rtype dict, None
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
        """Method to get a dictionary with the centroids of each cluster"""

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
        """Method to create the a dictionary with the ensembles (clustering with replacement)

        An rmsd and no. aligned residues threshold can be used, but if a qscore threshold is given, this one will
        be used instead

        :param rmsd_threshold: the rmsd threshold at which fragments are included into the ensemble (default 0.7)
        :type rmsd_threshold: float
        :param nalign_threshold: threshold no. aligned residues for a fragment to be included into an ensemble (default 30)
        :type nalign_threshold: int
        :param qscore_threshold: he qscore threshold at which fragments are included into the ensemble (default None)
        :type qscore_threshold: None, float
        :param nthreads: number of threads to compute the ensemble optimal parameter (default 1)
        :type nthreads: int
        :returns nothing
        :rtype None
        """

        self.logger.debug('Creating ensemble dictionary')
        self.ensemble_dict = {}

        # Check all centroids
        for clst_id in self.centroid_dict.keys():
            if clst_id != -1:
                self.ensemble_dict[clst_id] = self.get_ensemble(clst_id=clst_id, rmsd_threshold=rmsd_threshold,
                                                                nalign_threshold=nalign_threshold, nthreads=nthreads,
                                                                qscore_threshold=qscore_threshold)

    def get_ensemble(self, clst_id, rmsd_threshold=0.7, nalign_threshold=30, qscore_threshold=None, nthreads=1):
        """Search for fragments to form an ensemble given a cluster id

        Fragments within the threshold distance from the centroid will be included in the ensemble.
        Other centroids apart from the one corresponding with the cluster of interest are excluded from this search.

        :param clst_id: cluster id of the centroid of interest
        :type clst_id: int
        :param rmsd_threshold: the rmsd threshold at which fragments are included into the ensemble (default 0.7)
        :type rmsd_threshold: float
        :param nalign_threshold: threshold no. aligned residues for a fragment to be included into an ensemble (default 30)
        :type nalign_threshold: int
        :param qscore_threshold: he qscore threshold at which fragments are included into the ensemble (default None)
        :type qscore_threshold: None, float
        :param nthreads: number of threads to compute the ensemble optimal parameter (default 1)
        :type nthreads: int
        :returns hierarchy: pdb hierarchy with the ensemble (models are combined into an optimal structural alignment)
        :rtype :obj:`gemmi.Structure`
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

        self.logger.debug('Retrieving ensemble for cluster %s' % clst_id)
        fname = os.path.join(swamp.FRAG_PDB_DB, '{}.pdb')
        centroid = self.centroid_dict[clst_id]

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
            self.composition_dict[clst_id] = sorted(cluster_fragments)

        # If there are enough good fragments, make an ensemble and write into a file
        self.logger.debug('Retrieving optimum alignment for %s fragments [%s]' % (
            len(cluster_fragments), ', '.join(cluster_fragments)))
        gesamt, hierarchy = Gesamt.get_optimum_alignment([fname.format(frag_id) for frag_id in cluster_fragments],
                                                         nthreads=nthreads, logger=self.logger)
        if gesamt:
            self.logger.debug('Optimal qscore: %s' % gesamt.qscore)
        return hierarchy


class ParameterSearchResults(object):
    """Class to hold the results from a multi-threaded hyperparameter grid search

    Implements semaphore methods to regulate thread I/O into the result list

    :param :obj:`swamp.logger.swamplogger` logger: logger instance to record log messages
    :ivar :obj:`threading.lock` lock: lock to control I/O to the result instance
    :ivar :obj:`pandas.DataFrame` value: dataframe with the results of the grid search
    """

    def __init__(self, logger):
        self.lock = threading.Lock()
        self.value = pd.DataFrame(
            columns=["params_idx", "n_clst", "clst_size", "clst_qscore", "silhouette_score", "singletons"])
        self.logger = logger

    def register(self, new_results):
        self.logger.debug('Waiting for lock')
        self.lock.acquire()
        self.logger.debug('Acquired lock')
        self.value = pd.concat([self.value, new_results], 0)
        self.lock.release()
        self.logger.debug('Released lock')
