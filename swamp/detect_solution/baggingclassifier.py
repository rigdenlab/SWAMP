from swamp.detect_solution.classifier import Classifier
from scipy.stats import randint
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import BaggingClassifier
import pandas as pd
import matplotlib.pyplot as plt


# Class to hold decision tree ADA boost classifier
class DecisionTreeADABoost(Classifier):
    """Decision tree ADA boost classifier `"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {'n_estimators': randint(1000, 100000),
                'base_estimator__criterion': ['gini', 'entropy'],
                'base_estimator__splitter': ['best', 'random'],
                'base_estimator__class_weight': ['balanced', None],
                'base_estimator__max_features': randint(2, len(self.data.X_train.columns)),
                'base_estimator__min_samples_split': randint(2, 20),
                'base_estimator__max_depth': randint(5, 10),
                'base_estimator__min_samples_leaf': randint(1, 20),
                'base_estimator__max_leaf_nodes': randint(10, 20)}

    # Classifier to be trained
    @property
    def _clf(self):
        """Get the best hyperparams by random search`"""
        # create the decision tree
        clf1 = DecisionTreeClassifier(random_state=89)
        # create ada boost classifier
        return BaggingClassifier(base_estimator=clf1, random_state=22, bootstrap=True)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # get params for decission tree
        decission_tree_hyperparams = {}
        bagging_params = {}
        for hyperparam in self.best_params.keys():
            if "base_estimator" in hyperparam:
                decission_tree_hyperparams[hyperparam.split("__")[1]] = self.best_params[hyperparam]
            else:
                bagging_params[hyperparam] = self.best_params[hyperparam]
        # Create decission tree
        clf2 = DecisionTreeClassifier(random_state=89, **decission_tree_hyperparams)
        # Create adaboost
        self.trained_classifier = BaggingClassifier(clf2, random_state=22, bootstrap=True, n_jobs=self.n_jobs,
                                                    **bagging_params)
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
        feature_list = []
        for tree in self.trained_classifier.estimators_:
            feature_importances_ls = tree.feature_importances_
            feature_list.append(feature_importances_ls)

        df = pd.DataFrame(feature_list, columns=self.data.X_train.columns)
        df_mean = df[self.data.X_train.columns].mean(axis=0)
        df_mean.sort_values(ascending=False, inplace=True)
        df_std = df[self.data.X_train.columns].std(axis=0)
        df_mean.plot(kind='bar', colormap='viridis', yerr=[df_std], align="center", figsize=(20, 10), rot=90,
                     fontsize=16)
        plt.title('Histogram of Feature Importances for AdaBoost Classifier')
        plt.xlabel('Features')
        plt.tight_layout()
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()


if __name__ == '__main__':
    # Create the tree
    my_pred = DecisionTreeADABoost("/home/filo/Documents/Data/MR_edge_cases/ml_prepared_data_NEWGEN_balanced.csv")
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