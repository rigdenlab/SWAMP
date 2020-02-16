from swamp.clustering.clustering import Clustering
from sklearn.cluster import AffinityPropagation
from scipy.stats import randint, expon


class SwampAffinityPropagation(Clustering):
    """This class implements methods and datastructures to work with :py:obj:`sklearn.cluster.AffinityPropagation`

    :example:

    >>> import joblib
    >>> from swamp.clustering import SwampAffinityPropagation
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = SwampAffinityPropagation(dist_mtx)
    >>> my_clst.grid_search()
    >>> my_clst.cluster()
    >>> my_clst.assess_clustering(my_clst.labels)
    >>> my_clst.assess_clustering(my_clst.labels)
    """

    @property
    def _algorithm_name(self):
        """Name of the clustering algorithm (affinity)"""
        return "affinity"

    @property
    def _hyper_params(self):
        """Dictionary with the range of possible values for each of the clustering hyper-parameters"""

        return {"damping": expon(0.1),
                "convergence_iter": randint[15, 100]
                }

    def _clustering(self, **kwargs):
        """Perform clustering with a given set of arguments"""
        return AffinityPropagation(affinity='precomputed', **kwargs)

    def cluster(self):
        """Method to perform a clustering using the :py:attr:`~swamp.clustering.Clustering.best_params`

        :raises ValueError: the attribute :py:attr:`~swamp.clustering.Clustering.similarity_mtx` is None
        """

        self.logger.info(self.clustering_header)

        if self.similarity_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = AffinityPropagation(affinity='precomputed', **self.best_params)
        clst.fit(self.similarity_mtx)
        self.labels = clst.labels_
