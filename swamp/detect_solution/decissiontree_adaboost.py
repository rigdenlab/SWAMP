from swamp.detect_solution.classifier import Classifier
from scipy.stats import randint
from scipy.stats import uniform
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
import pandas as pd
import matplotlib.pyplot as plt


# Class to hold decision tree ADA boost classifier
class DecisionTreeADABoost(Classifier):
    """Decision tree ADA boost classifier `"""

    # set up hyperparameters
    @property
    def _hyper_params(self):
        return {'n_estimators': randint(1000, 100000),
                'learning_rate': uniform(0.001, 1.0),
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
        clf1 = DecisionTreeClassifier(random_state=65)
        # create ada boost classifier
        return AdaBoostClassifier(base_estimator=clf1, random_state=100)

    # Train the classifier
    def train_classifier(self, get_stats=False):
        # get params for decission tree
        decission_tree_hyperparams = {}
        ada_params = {}
        for hyperparam in self.best_params.keys():
            if "base_estimator" in hyperparam:
                decission_tree_hyperparams[hyperparam.split("__")[1]] = self.best_params[hyperparam]
            else:
                ada_params[hyperparam] = self.best_params[hyperparam]
        # Create decission tree
        clf2 = DecisionTreeClassifier(random_state=0, **decission_tree_hyperparams)
        # Create adaboost
        self.trained_classifier = AdaBoostClassifier(clf2, **ada_params)
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
    my_pred = DecisionTreeADABoost(
        "/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_prepared_data_NEWGEN_unbalanced.csv")
    my_pred.n_jobs = 64
    my_pred.n_iter = 1000
    my_pred.cv_folds = 5
    print(my_pred.data.dataframe.shape)

    # Balance data
    my_pred.balance_data()

    # Omit some features
    my_pred.omit_samples(feature="MODIF", omit_value="singlechain")
    my_pred.omit_features(["TARGET_ID", "CLST_ID", "CON_SCO", "N_MODELS", "MODIF", "INIT_RFACT", "FINAL_RFACT"])

    # REFMAC FEATURES
    my_pred.encode_deltas(initial_field="INIT_RFREE", final_field="FINAL_RFREE", new_fieldname="RFREE_DELTA")
    my_pred.encode_deltas(initial_field="INIT_CHIRVOL", final_field="FINAL_CHIRVOL", new_fieldname="CHIRVOL_DELTA")
    my_pred.encode_deltas(initial_field="INIT_BONDLENGHT", final_field="FINAL_BONDLENGHT",
                          new_fieldname="BONDLENGHT_DELTA")
    my_pred.encode_deltas(initial_field="INIT_BONDANGLE", final_field="FINAL_BONDANGLE",
                          new_fieldname="BONDANGLE_DELTA")

    # SPACEGROUP POLARITY
    my_pred.encode_spacegroup(insitu=True)

    # Split dataset
    my_pred.split_dataset()
    print(my_pred.features)
    print(my_pred.data.dataframe.shape)
    # Random search
    my_pred.random_hyperparameter_search()
    # Training
    my_pred.train_classifier(get_stats=True)
    my_pred.feature_importance()
    # Prediction
    my_pred.predict(X=my_pred.data.X_test, get_stats=True)
    # Plots
    my_pred.draw_conf_matrix(y_test=my_pred.data.y_test, y_pred=my_pred.y_pred)
    my_pred.plot_precision_recall_vs_threshold(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
    my_pred.plot_roc_curve(y_test=my_pred.data.y_test, y_test_class_one=my_pred.y_pred_proba[:, 1])
