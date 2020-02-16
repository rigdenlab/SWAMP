from swamp.clustering.clustering import Clustering
from sklearn.cluster import SpectralClustering
from scipy.stats import randint, expon


class Spectral(Clustering):
    """This class implements methods and datastructures to work with :py:obj:`sklearn.cluster.Spectral`

    :example:

    >>> from swamp.clustering import Spectral
    >>> import joblib
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = Spectral(dist_mtx)
    >>> my_clst.grid_search()


    """

    @property
    def _algorithm_name(self):
        """Name of the clustering algorithm (spectral)"""
        return "spectral"

    @property
    def _hyper_params(self):
        """Dictionary with the range of possible values for each of the clustering hyper-parameters"""

        return {"n_clusters": randint(200, 900),
                "eigen_solver": [None, "arpack", "lobpcg"],
                "assign_labels": ["kmeans", "discretize"],
                "n_neighbors": randint(2, 10),
                "gamma": expon(0.1),
                }

    def _clustering(self, **kwargs):
        """Perform clustering with a given set of arguments"""
        return SpectralClustering(affinity='precomputed', n_jobs=1, **kwargs)

    def cluster(self):
        """Method to perform a clustering using the :py:attr:`~swamp.clustering.Clustering.best_params`

        :raises ValueError: the attribute :py:attr:`~swamp.clustering.Clustering.similarity_mtx` is None
        """

        self.logger.info(self.clustering_header)

        if self.similarity_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = SpectralClustering(n_jobs=self.nthreads, affinity='precomputed', **self.best_params)
        clst.fit(self.similarity_mtx)
        self.labels = clst.labels_
