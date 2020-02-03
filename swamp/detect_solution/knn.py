from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from swamp.detect_solution.classifier import Classifier
from scipy.stats import randint
from sklearn.model_selection import cross_val_score
import numpy as np


# Class to hold decision tree ADA boost classifier
class KNN(Classifier):
    """KNNs classifier `"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {'n_neighbors': randint(5, 100),
                'weights': ["uniform", "distance"],
                'algorithm': ["auto", "ball_tree", "kd_tree", "brute"],
                'leaf_size': randint(10, 200),
                'p': [1, 2]
                }

    # Classifier to be trained
    @property
    def _clf(self):
        return KNeighborsClassifier()

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # Scale data
        if not self.data.is_scaled:
            self.data._scale()
        # Create KNN
        self.trained_classifier = KNeighborsClassifier(n_jobs=self.n_jobs, **self._best_params)
        self.trained_classifier.fit(self.data.X_train, self.data.y_train)
        # Get stats
        if get_stats:
            self.get_training_stats()
        # now this is trained
        self.is_trained = True

    def predict(self, X, get_stats=False):
        """Predict the class of a set of given features"""
        self.y_pred = self.trained_classifier.predict(X)
        self.y_pred_proba = self.trained_classifier.predict_proba(X)
        # get stats
        if get_stats:
            self.get_prediction_stats(y_test=self.data.y_test, y_pred=self.y_pred)

    def feature_importance(self, fname=None):
        """Get the feature importance after the classification is made`"""
        # train the classifier
        if not self.is_trained:
            self.train_classifier()
        feature_score_list = []
        for i in range(len(self.features)):
            X = self.data.X_train_unscaled.iloc[:, i].values.reshape(-1, 1)
            scores = cross_val_score(self.trained_classifier, X, self.data.y_train, cv=self.cv_folds)
            feature_score_list.append(scores.mean())
        feature_importance = sorted(zip(feature_score_list, self.features), reverse=True)
        print("Feature importance:\n{}".format(feature_importance))
        feature = list(zip(*feature_importance))[1]
        score = list(zip(*feature_importance))[0]
        x_pos = np.arange(len(feature))
        plt.figure()
        plt.bar(x_pos, score, align='center', color='b')
        plt.xticks(x_pos, feature, rotation=90, fontsize=16)
        plt.title('Histogram of Feature Importances for best KNN')
        plt.xlabel('Features')
        plt.tight_layout()
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()

if __name__ == '__main__':
    my_pred = KNN("/home/filo/Documents/Data/MR_edge_cases/ml_prepared_data_NEWGEN_balanced.csv")
    my_pred.n_jobs = 12
    my_pred.n_iter = 1
    my_pred.omit_samples(feature="TARGET_ID", omit_value="3zot")
    my_pred.omit_features(["TARGET_ID", "CLST_ID", "CON_SCO", "MODIF", "eLLG", "N_MODELS"])
    my_pred.balance_data()
    my_pred.split_dataset()
    my_pred.random_hyperparameter_search()
    my_pred.train_classifier()
    my_pred.feature_importance()
    my_pred.predict(X=my_pred.data.X_test, get_stats=True)
    my_pred.draw_conf_matrix(y_test=my_pred.data.y_test, y_pred=my_pred.y_pred)
    my_pred.plot_precision_recall_vs_threshold(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
    my_pred.plot_roc_curve(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
