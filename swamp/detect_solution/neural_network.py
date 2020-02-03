from sklearn.neural_network import MLPClassifier
from swamp.detect_solution.classifier import Classifier
from scipy.stats import expon, randint, uniform


# Class to hold decision tree ADA boost classifier
class NeuralNetwork(Classifier):
    """Multi-layer Perceptron classifier"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {'hidden_layer_sizes': randint(50, 300),
                'activation': ['identity', 'logistic', 'tanh', 'relu'],
                'solver': ["lbfgs", "sgd", "adam"],
                'alpha': expon(scale=0.01),
                'learning_rate': ['constant', 'invscaling', 'adaptative'],
                'learning_rate_init': expon(0.1),
                'power_t': expon(scale=1.0),
                'max_iter': randint(200, 1000),
                'momentum': uniform(0, 1),
                'nesterovs_momentum': [True, False],
                'beta_1': uniform(0, 1),
                'beta_2': uniform(0, 1),
                'epsilon': expon(0.00001),
                }

    # Classifier to be trained
    @property
    def _clf(self):
        # create the svm
        return MLPClassifier(random_state=57)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # Scale data
        if not self.data.is_scaled:
            self.data._scale()
        # Create SVM
        self.trained_classifier = MLPClassifier(random_state=57, **self._best_params)
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
        print("Feature importance is not available for Neural networks!")


if __name__ == '__main__':
    my_pred = NeuralNetwork("/home/filo/Documents/Data/MR_edge_cases/ml_prepared_data_NEWGEN_balanced.csv")
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
