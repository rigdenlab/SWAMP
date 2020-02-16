from swamp.clustering.clustering import Clustering
from sklearn.cluster import OPTICS
from scipy.stats import randint


class Optics(Clustering):
    """This class implements methods and datastructures to work with :py:obj:`sklearn.cluster.OPTICS`

    :example:

    >>> import joblib
    >>> from swamp.clustering import Optics
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = Optics(dist_mtx)
    >>> my_clst.grid_search()
    >>> my_clst.cluster()
    >>> my_clst.assess_clustering(my_clst.labels)
    >>> my_clst.assess_clustering(my_clst.labels)
    """

    @property
    def _algorithm_name(self):
        """Name of the clustering algorithm (optics)"""
        return "optics"

    @property
    def _hyper_params(self):
        """Dictionary with the range of possible values for each of the clustering hyper-parameters"""

        return {"min_samples": randint(2, 5),
                "max_eps": [0.5, 0.6, 0.7, 0.8, 0.9],
                "p": [1, 2],
                "cluster_method": ['xi', 'dbscan'],
                "xi": [0.05, 0.01, 0.1, 0.02, 0.03, 0.07],
                "predecessor_correction": [True, False],
                "min_cluster_size": [2, 3],
                "leaf_size": randint(10, 50),
                }

    def _clustering(self, **kwargs):
        """Perform clustering with a given set of arguments"""
        return OPTICS(metric='precomputed', n_jobs=1, **kwargs)

    def cluster(self):
        """Method to perform a clustering using the :py:attr:`~swamp.clustering.Clustering.best_params`

        :raises ValueError: the attribute :py:attr:`~swamp.clustering.Clustering.similarity_mtx` is None
        """

        self.logger.info(self.clustering_header)

        if self.distance_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = OPTICS(metric='precomputed', n_jobs=self.nthreads, **self.best_params)
        clst.fit(self.distance_mtx)
        self.labels = clst.labels_
