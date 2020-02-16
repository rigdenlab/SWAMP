from swamp.clustering import Clustering
from sklearn.cluster import SpectralClustering
from scipy.stats import randint, expon


class Spectral(Clustering):
    """ Class to wrap a sklearn spectral clustering instance

    :example

    >>> from swamp.clustering import Spectral
    >>> import joblib
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = Spectral(dist_mtx)
    >>> my_clst.parameter_search()


    """

    @property
    def _algorithm_name(self):
        return "spectral"

    @property
    def _hyper_params(self):
        """Property containing all the hyper parameters that can be modified for spectral clustering"""

        return {"n_clusters": randint(200, 900),
                "eigen_solver": [None, "arpack", "lobpcg"],
                "assign_labels": ["kmeans", "discretize"],
                "n_neighbors": randint(2, 10),
                "gamma": expon(0.1),
                }

    def _clustering(self, **kwargs):
        return SpectralClustering(affinity='precomputed', n_jobs=1, **kwargs)

    def cluster(self):
        """Method to perform a clustering using the best parameters and over a given distance matrix"""

        self.logger.info(self.clustering_header)

        if self.similarity_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = SpectralClustering(n_jobs=self.nthreads, affinity='precomputed', **self.best_params)
        clst.fit(self.similarity_mtx)
        self.labels = clst.labels_
