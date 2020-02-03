from sklearn.svm import SVC
import matplotlib.pyplot as plt
from swamp.detect_solution.classifier import Classifier
from scipy.stats import expon
import pandas as pd


# Class to hold decision tree ADA boost classifier
class SVM(Classifier):
    """SVM classifier `"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {'C': expon(100),
                'shrinking': [True, False],
                'decision_function_shape': ["ovo", "ovr"],
                'class_weight': ['balanced', None]
                }

    # Classifier to be trained
    @property
    def _clf(self):
        # create the svm
        return SVC(kernel="linear", probability=True, random_state=100)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # Scale data
        if not self.data.is_scaled:
            self.data._scale()
        # Create SVM
        self.trained_classifier = SVC(kernel="linear", probability=True, random_state=100, **self._best_params)
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
        feature_importances_df = pd.DataFrame()
        feature_importances_df['importance'] = [x for x in self.trained_classifier.coef_.tolist()[0]]
        feature_importances_df['feature'] = [x for x in self.data.X_train_unscaled.columns]
        feature_importances_df.sort_values(by='importance', ascending=False, inplace=True)
        feature_importances_df = feature_importances_df.set_index('feature', drop=True)
        plt.figure()
        feature_importances_df.plot.bar(figsize=(20, 10), colormap='viridis', fontsize=16)
        plt.title('Histogram of Feature Importances for SVM Classifier')
        plt.xlabel('Features')
        plt.tight_layout()
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()


if __name__ == '__main__':
    my_pred = SVM("/home/filo/Documents/Data/MR_edge_cases/ml_prepared_data_NEWGEN_balanced.csv")
    my_pred.n_jobs = 12
    my_pred.n_iter = 1
    my_pred.omit_samples(feature="TARGET_ID", omit_value="3zot")
    my_pred.balance_data()
    my_pred.omit_features(["TARGET_ID", "CLST_ID", "CON_SCO", "MODIF", "eLLG", "N_MODELS"])
    my_pred.split_dataset()
    my_pred.random_hyperparameter_search()
    my_pred.train_classifier()
    my_pred.feature_importance()
    my_pred.predict(X=my_pred.data.X_test, get_stats=True)
    my_pred.draw_conf_matrix(y_test=my_pred.data.y_test, y_pred=my_pred.y_pred)
    my_pred.plot_precision_recall_vs_threshold(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
    my_pred.plot_roc_curve(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
