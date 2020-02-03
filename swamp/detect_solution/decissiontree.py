from sklearn.tree import DecisionTreeClassifier
import matplotlib.pyplot as plt
from swamp.detect_solution.classifier import Classifier
from scipy.stats import randint, expon
import pandas as pd


# Class to hold decision tree ADA boost classifier
class RandomForests(Classifier):
    """Random Forests classifier `"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {"criterion": ["gini", "entropy"],
                'splitter':["best", "random"],
                'max_features': randint(2, len(self.data.X_train.columns)), #max number of features when splitting
                "min_samples_split": randint(2, 20),#min samples per node to induce split
                "max_depth": randint(3, 10),#max number of splits to do
                "min_samples_leaf": randint(1, 20),#min number of samples in a leaf; may set to 1 anyway
                "max_leaf_nodes": randint(10, 20) #max number of leaves
                }

    # Classifier to be trained
    @property
    def _clf(self):
        return DecisionTreeClassifier(random_state=36)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # Create Random forest
        self.trained_classifier = DecisionTreeClassifier(random_state=36, **self._best_params)
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

        feature_importances_df = pd.DataFrame({'importance': self.trained_classifier.feature_importances_})
        feature_importances_df['feature'] = self.data.X_train.columns
        feature_importances_df.sort_values(by='importance', ascending=False, inplace=True)
        feature_importances_df = feature_importances_df.set_index('feature', drop=True)
        plt.figure()
        feature_importances_df.plot.bar(figsize=(20, 10), colormap='viridis', fontsize=16)
        plt.title('Histogram of Feature Importances for Decission tree Classifier')
        plt.xlabel('Features')
        plt.tight_layout()
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()


if __name__ == '__main__':
    my_pred = RandomForests("/home/filo/Documents/Data/MR_edge_cases/ml_prepared_data_NEWGEN_balanced.csv")
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
