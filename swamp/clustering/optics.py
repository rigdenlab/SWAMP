from swamp.clustering import Clustering
from sklearn.cluster import OPTICS
from scipy.stats import randint


class Optics(Clustering):
    """ Class to wrap a sklearn Optics clustering instance

    :example

    >>> import joblib
    >>> from swamp.clustering import Optics
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = Optics(dist_mtx)
    >>> my_clst.parameter_search()
    >>> my_clst.cluster()
    >>> my_clst.assess_clustering(my_clst.labels)
    >>> my_clst.assess_clustering(my_clst.labels)
    """

    @property
    def _algorithm_name(self):
        return "optics"

    @property
    def _hyper_params(self):
        """Property containing all the hyper parameters that can be modified for optics clustering"""

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
        return OPTICS(metric='precomputed', n_jobs=1, **kwargs)

    def cluster(self):
        """Method to perform a clustering using the best parameters and over a given distance matrix"""

        self.logger.info(self.clustering_header)

        if self.distance_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = OPTICS(metric='precomputed', n_jobs=self.nthreads, **self.best_params)
        clst.fit(self.distance_mtx)
        self.labels = clst.labels_
