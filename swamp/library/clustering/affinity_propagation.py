from swamp.library.clustering.clustering import Clustering
from sklearn.cluster import AffinityPropagation
from scipy.stats import randint, expon


class AffProp(Clustering):
    """ Class to wrap a sklearn AffinityPropagation clustering instance

    :example

    >>> import joblib
    >>> from swamp.library.clustering.affinity_propagation import AffProp
    >>> dist_mtx = joblib.load('<dist_mtx.pckl>')
    >>> dist_mtx = dist_mtx.fillna(0)
    >>> my_clst = AffProp(dist_mtx)
    >>> my_clst.parameter_search()
    >>> my_clst.cluster()
    >>> my_clst.assess_clustering(my_clst.labels)
    >>> my_clst.assess_clustering(my_clst.labels)
    """

    @property
    def _algorithm_name(self):
        return "affinity"

    @property
    def _hyper_params(self):
        """Property containing all the hyper parameters that can be modified for affinity propagation clustering"""

        return {"damping": expon(0.1),
                "convergence_iter": randint[15, 100]
                }

    def _clustering(self, **kwargs):
        return AffinityPropagation(affinity='precomputed', **kwargs)

    def cluster(self):
        """Method to perform a clustering using the best parameters and over a given distance matrix"""

        self.logger.info(self.clustering_header)

        if self.similarity_mtx is None:
            raise ValueError('Need to load a distance matrix before clustering!')

        clst = AffinityPropagation(affinity='precomputed', **self.best_params)
        clst.fit(self.similarity_mtx)
        self.labels = clst.labels_
