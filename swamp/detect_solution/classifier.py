import abc
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import cross_val_score
from sklearn import metrics
from swamp.detect_solution.data import Data
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

ABC = abc.ABCMeta('ABC', (object,), {})


# Create a class to contain the classifier
class Classifier(ABC):
    # Constructor method
    def __init__(self, csv_data):
        self.data = Data(csv_data)
        self.y_pred = None
        self.y_pred_proba = None
        self._n_jobs = 1
        self._n_iter = 10
        self._cv_folds = 3
        self._scoring_metric = "accuracy"
        self._best_params = None
        self._is_trained = False
        self.trained_classifier = None

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def train_classifier(self, get_stats=False):
        """ Abstract method to train a given Classifier"""
        pass

    @abc.abstractmethod
    def predict(self, X, get_stats=False):
        """Abstract method to make a prediction using Classifier"""
        pass

    @abc.abstractmethod
    def feature_importance(self, fname):
        """Abstract method to make a prediction using Classifier"""
        pass

    @property
    @abc.abstractmethod
    def _hyper_params(self):
        """Abstract property to contain hyperparams to optimise"""
        pass

    @property
    @abc.abstractmethod
    def _clf(self):
        """Abstract property to contain hyperparams to optimise"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def features(self):
        """Property with the colnames in X"""
        if self.data.X is None:
            return None
        else:
            return self.data.X_train_unscaled.columns

    @property
    def n_jobs(self):  # pragma: no cover
        """Property with the number of available parallel jobs for execution of :obj:`~pyjob.task.Task`"""
        return self._n_jobs

    @n_jobs.setter
    def n_jobs(self, value):
        self._n_jobs = value

    @property
    def is_trained(self):  # pragma: no cover
        """Property to determine if the classifier is already trained or not :obj:`~pyjob.task.Task`"""
        return self._is_trained

    @is_trained.setter
    def is_trained(self, value):
        self._is_trained = value

    @property
    def n_iter(self):  # pragma: no cover
        """Property with the number of iterations for the randomsearch of for execution of :obj:`~pyjob.task.Task`"""
        return self._n_iter

    @n_iter.setter
    def n_iter(self, value):
        self._n_iter = value

    @property
    def cv_folds(self):  # pragma: no cover
        """Property with the number of iterations for the randomsearch of for execution of :obj:`~pyjob.task.Task`"""
        return self._cv_folds

    @cv_folds.setter
    def cv_folds(self, value):
        self._cv_folds = value

    @property
    def scoring_metric(self):  # pragma: no cover
        """Property with the scoring metrix for the random search of :obj:`~pyjob.task.Task`"""
        return self._scoring_metric

    @scoring_metric.setter
    def scoring_metric(self, value):
        self._scoring_metric = value

    @property
    def best_params(self):  # pragma: no cover
        """Property with the best params for execution of :obj:`~pyjob.task.Task`"""
        return self._best_params

    @best_params.setter
    def best_params(self, value):
        self._best_params = value

    @property
    def _info_shape(self):
        return self.data._info_shape

    # ------------------ Some general methods ------------------

    # Create a method to save the classifer into a pickle
    def save_pickle(self, fname):
        joblib.dump(self.trained_classifier, fname)

    # Create a method to load a trained classifier from a pickle
    def load_pickle(self, fname):
        self.trained_classifier = joblib.load(fname)

    # Create a method to ommit samples with NA in a given list of features
    def omit_NAs(self, feature_list):
        self.data._omit_NAs(feature_list)

    # Create a method to ommit a list of features from the dataframe
    def omit_features(self, feature_list):
        self.data._omit_features(feature_list)

    # Create a method to omit samples matching a patters (e.g. target_id)
    def omit_samples(self, feature, omit_value):
        self.data._omit_samples(feature, omit_value)

    # Create a function for one hot encoding of a particular feature
    def one_hot_encode(self, feature):
        self.data._one_hot_encode(feature)

    # Function to encode the spacegroup info
    def encode_spacegroup(self, **kwargs):
        self.data._encode_spacegroup(**kwargs)

    # Method to get polynomial interactions
    def create_polynomials(self, **kwargs):
        self.data._create_polynomials(**kwargs)

    # Method to get training and testing datasets
    def split_dataset(self, **kwargs):
        self.data._split_dataset(**kwargs)

    # Method to get the calibration set seprated
    def get_calibration_dataset(self, **kwargs):
        self.data._extract_calibration_set(**kwargs)

    # Function to scale X
    def scale_data(self):
        self.data._scale()

    # Function to balance the classes in the data
    def balance_data(self):
        self.data._balance_dataset()

    # Function to encode the deltas
    def encode_deltas(self, **kwargs):
        self.data._encode_delta(**kwargs)

    # Perform a randomised search for optimal hyperparameters
    def random_hyperparameter_search(self):
        """Method to make a random search of hyperparams for Classifier"""
        # Make sure the user split the data
        if not self.data.is_split:
            self.split_dataset()
        # Maybe it's necessary to scale data
        if isinstance(self._clf, KNeighborsClassifier) or isinstance(self._clf, SVC) and not self.data.is_scaled:
            self.scale_data()
        # Do the random search
        rand_search = RandomizedSearchCV(self._clf, self._hyper_params, random_state=5, cv=self.cv_folds,
                                         n_iter=self.n_iter, scoring=self.scoring_metric, n_jobs=self.n_jobs)
        rand_search = rand_search.fit(self.data.X_train, self.data.y_train)
        print('Best parameters: ' + str(rand_search.best_params_) + '\n')
        print('Best score: ' + str(rand_search.best_score_) + '\n')
        # save best params
        self.best_params = rand_search.best_params_

    def grid_hyperparameter_search(self):
        """Method to make a grid search of hyperparams for the Classifier"""
        # build a grid searcher
        grid_search = GridSearchCV(self._clf, self._hyper_params, scoring=self.scoring_metric, n_jobs=self.n_jobs, cv=3)
        grid_search.fit(self.data.X_train, self.data.y_train)
        print('Best parameters: ' + str(grid_search.best_params_) + '\n')
        print('Best score: ' + str(grid_search.best_score_) + '\n')
        # save best params
        self.best_params = grid_search.best_params_

    # Create a method to get some basic stats from the trained classifier
    def get_training_stats(self):
        """Method to get the stats for Classifier`"""
        # distribution --> accuracy
        accuracy_each_cv = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                           scoring='accuracy')
        accuracy_mean_cv = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                           scoring='accuracy').mean()
        # calculate cross_val_scoring with different scoring functions for CV train set
        train_roc_auc = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                        scoring='roc_auc').mean()
        train_accuracy = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                         scoring='accuracy').mean()
        train_recall = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                       scoring='recall').mean()
        train_precision = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                          scoring='precision').mean()
        train_f1 = cross_val_score(self.trained_classifier, self.data.X_train, self.data.y_train, cv=3,
                                   scoring='f1').mean()

        print('Get various cross_val_scores to evaluate clf performance for best parameters \n')
        print('Accuracy for each of 3 CV folds: %s \n' % accuracy_each_cv)
        print('Mean accuracy over all 3 CV folds: %s \n' % accuracy_mean_cv)
        print('ROC_AUC mean for 3-fold CV: %s \n' % train_roc_auc)
        print('Accuracy mean for 3-fold CV: %s \n' % train_accuracy)
        print('Recall mean for 3-fold CV: %s \n' % train_recall)
        print('Precision mean for 3-fold CV: %s \n' % train_precision)
        print('F1 score mean for 3-fold CV: %s \n' % train_f1)

    # Method to get test prediction stats
    @staticmethod
    def get_prediction_stats(y_test, y_pred):
        # calculate accuracy
        y_accuracy = metrics.accuracy_score(y_test, y_pred)
        # examine the class distribution of the testing set (using a Pandas Series method)
        class_dist = y_test.value_counts()
        # calculate the percentage of ones
        # because y_test only contains ones and zeros, we can simply calculate the mean = percentage of ones
        ones = y_test.mean()
        # calculate the percentage of zeros
        zeros = 1 - y_test.mean()
        # calculate null accuracy in a single line of code
        # only for binary classification problems coded as 0/1
        null_acc = max(y_test.mean(), 1 - y_test.mean())
        print('Accuracy score or agreement between y_test and y_pred_class: %s \n' % y_accuracy)
        print('Class distribution for y_test:\n%s\n' % class_dist)
        print('Percent 1s in y_test: %s \n' % ones)
        print('Percent 0s in y_test: %s \n' % zeros)
        print('Null accuracy in y_test: %s \n' % null_acc)

    # Method to draw confussion matrix
    @staticmethod
    def draw_conf_matrix(y_test, y_pred, fname=None):
        plt.figure()
        conf_mat_test = metrics.confusion_matrix(y_test, y_pred)
        labels = ['0', '1']
        ax = plt.subplot()
        sns.heatmap(conf_mat_test, annot=True, ax=ax, fmt="d", cmap="viridis")
        plt.title('Confusion matrix of the classifier')
        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)
        plt.xlabel('Predicted')
        plt.ylabel('True')
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()

    # plot precision and recall curve
    @staticmethod
    def plot_precision_recall_vs_threshold(y_test, y_test_class_one, fname=None):
        """ y_test_class_one: test data to be class 1 (  y_test_class_one = y_pred_proba[:, 1]  )"""
        plt.figure()
        precisions, recalls, thresholds_tree = metrics.precision_recall_curve(y_test, y_test_class_one)
        plt.plot(thresholds_tree, precisions[:-1], "b--", label="Precision")
        plt.plot(thresholds_tree, recalls[:-1], "g--", label="Recall")
        plt.title('Precsion-Recall plot for for fragment swamp success classifier using class 1')
        plt.xlabel("Threshold")
        plt.legend(loc="upper left")
        plt.ylim([0, 1])
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()

    # plot ROC
    @staticmethod
    def plot_roc_curve(y_test, y_test_class_one, fname=None):
        plt.figure()
        AUC_test_class1 = metrics.roc_auc_score(y_test, y_test_class_one)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_test_class_one)
        plt.plot(fpr, tpr, linewidth=2)
        plt.plot([0, 1], [0, 1], 'k--')
        plt.axis([0, 1, 0, 1])
        plt.title('ROC curve for fragment swamp success classifier for class 1 (AUC = %s)' % round(AUC_test_class1, 2))
        plt.xlabel('False Positive Rate (1 - Specificity)')
        plt.ylabel('True Positive Rate (Sensitivity)')
        plt.grid(True)
        if fname is not None:
            plt.savefig(fname)
        else:
            plt.show()
