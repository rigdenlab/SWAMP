from sklearn.ensemble import RandomForestClassifier
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
        return {'n_estimators': randint(100, 10000),
                'criterion': ["gini", "entropy"],
                'max_depth': randint(5, 10),
                'min_samples_split': randint(2, 20),
                'min_samples_leaf': randint(1, 20),
                #'min_weight_fraction_leaf': randint(0, 0.5),
                'max_features': randint(2, len(self.data.X_train.columns)),
                'max_leaf_nodes': randint(10, 200),
                'min_impurity_decrease': expon(0.1),
                'class_weight': ['balanced', None]
                }

    # Classifier to be trained
    @property
    def _clf(self):
        return RandomForestClassifier(random_state=89)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # Create Random forest
        self.trained_classifier = RandomForestClassifier(random_state=89, n_jobs=self.n_jobs, **self._best_params)
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
        feature_list = []
        for tree in self.trained_classifier.estimators_:
            feature_importances_ls = tree.feature_importances_
            feature_list.append(feature_importances_ls)
        df = pd.DataFrame(feature_list, columns=self.features)
        df_mean = df[self.features].mean(axis=0)
        df_mean.sort_values(ascending=False, inplace=True)
        df_std = df[self.features].std(axis=0)
        plt.figure()
        df_mean.plot(kind='bar', colormap='viridis', yerr=[df_std], align="center", figsize=(20, 10), rot=90, fontsize=16)
        plt.title('Histogram of Feature Importances for Random Forest Classifier')
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
