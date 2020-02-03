import joblib
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PolynomialFeatures

_POLAR_SPs = ("P1", "P2", "P21", "C2", "P4", "P41", "P42", "P43", "I4", "I41", "P3", "P31", "P32", "R3", "P6", "P61",
              "P65", "P62", "P64", "P63", "P1211", "C121")
_MULTIPLE_INDEXING_SPs = ("P312", "P321", "P3112", "P3121", "P23", "F23", "I23", "P213", "I213", "F23")


# Create a class to hold the data contained in a given csv
class Data(object):

    # Constructor method
    def __init__(self, csv_file):
        self.fname = csv_file
        self.dataframe = pd.read_csv(csv_file, na_filter=False, skipinitialspace=True, thousands=',')
        self.X = None
        self.y = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.is_scaled = False
        self.X_train_unscaled = None
        self.is_split = False
        self.X_test_unscaled = None
        self.is_balanced = False
        self.X_calibration = None
        self.y_calibration = None
        self.X_calibration_unscaled = None

    @property
    def _info_shape(self):
        return {"SOL==0": self.dataframe[self.dataframe.SOLUTION == 0].shape,
                "SOL==1": self.dataframe[self.dataframe.SOLUTION == 1].shape,
                "ALL": self.dataframe.shape}

    # Create a method to save the data into a pickle
    def save_pickle(self, fname):
        joblib.dump(self.dataframe, fname)

    # Create a method to load a data classifier from a pickle
    def load_pickle(self, fname):
        self.dataframe = joblib.load(fname)

    # Method to sort the dataframe by a given feature
    def _sort_dataframe(self, by_field, reset_index=True):
        self.dataframe = self.dataframe.sort_values(by=by_field, ascending=False)
        if reset_index:
            self.dataframe = self.dataframe.reset_index()

    # Method to catch rows mattching a given pattern in the dataframe
    def _catch_samples(self, field, target):
        self.dataframe = self.dataframe[self.dataframe[field] == target]

    # Create a method to ommit samples with NA in a given list of features
    def _omit_NAs(self, feature_list):
        for field in feature_list:
            self.dataframe = self.dataframe[self.dataframe[str(field)] != "NA"]

    # Create a function to substitute "NA" with np.nan
    def _replace_nan(self):
        self.dataframe = self.dataframe.replace("NA", np.nan)

    # Create a method to ommit a list of features from the dataframe
    def _omit_features(self, feature_list):
        for field in feature_list:
            self.dataframe = self.dataframe.drop(str(field), 1)

    # Create a method to omit samples matching a patters (e.g. target_id)
    def _omit_samples(self, feature, omit_value):
        self.dataframe = self.dataframe[self.dataframe[feature] != omit_value]

    # Create a function for one hot encoding of a particular feature
    def _one_hot_encode(self, feature):
        one_hot = pd.get_dummies(self.dataframe[str(feature)])
        self.dataframe = self.dataframe.join(one_hot)
        self.dataframe = self.dataframe.drop(str(feature), 1)

    # Function to encode the delta information
    def _encode_delta(self, initial_field, final_field, new_fieldname, percentage=False, drop=False):
        if not percentage:
            self.dataframe[new_fieldname] = self.dataframe[final_field].astype(float) - self.dataframe[
                initial_field].astype(float)
        else:
            self.dataframe["tmp_diff"] = self.dataframe[final_field] - self.dataframe[initial_field]
            self.dataframe[new_fieldname] = self.dataframe["tmp_diff"] / self.dataframe[initial_field]
            self.dataframe = self.dataframe.drop("tmp_diff", 1)
        if drop:
            self.dataframe = self.dataframe.drop(initial_field, 1)
            self.dataframe = self.dataframe.drop(final_field, 1)

    # Function to encode the spacegroup information
    def _encode_spacegroup(self, insitu=False, colname="SPACEGROUP", new_colname="SPACEGROUP_POLAR",
                           target_spacegroups=_POLAR_SPs):
        # Make a new field
        if not insitu:
            self.dataframe[new_colname] = self.dataframe[colname]
            colname = new_colname
        # Change values
        self.dataframe.loc[self.dataframe[colname].isin(target_spacegroups), [colname]] = 1
        self.dataframe.loc[self.dataframe[colname] != 1, [colname]] = 0

    # Function to create a new polynomial set of features [1, a, b, a^2, ab, b^2]
    def _create_polynomials(self, feature_list, drop=True, degree=2, interaction_only=False):
        # Get the features for the polynomials and get the polynomials
        feature_slice = self.dataframe[feature_list]
        poly = PolynomialFeatures(degree=degree, include_bias=False, interaction_only=interaction_only)
        polynomial_features = pd.DataFrame(poly.fit_transform(feature_slice))

        # We need to get the new colnames into a dict
        feature_template = "x{}"
        squared_feature_template = "{}^2"
        feature_dict = {}
        for idx, feature in enumerate(feature_list):
            feature_dict[feature_template.format(idx)] = feature
            feature_dict[
                squared_feature_template.format(feature_template.format(idx))] = squared_feature_template.format(
                feature)
        # Get the new colnames
        new_colnames = []
        new_feature_template = '{}_x_{}'
        for feature in poly.get_feature_names():
            if ' ' not in feature:
                new_colnames.append(feature_dict[feature])
            else:
                feature = feature.split()
                new_colnames.append(new_feature_template.format(feature_dict[feature[0]], feature_dict[feature[1]]))
        # Get rid of the original feature list
        for feature in feature_list:
            self.dataframe = self.dataframe.drop(feature, 1)
        polynomial_features.columns = new_colnames
        # Merge the DFs
        self.dataframe = self.dataframe.reset_index(drop=True)
        polynomial_features = polynomial_features.reset_index(drop=True)
        self.dataframe = pd.concat([self.dataframe, polynomial_features], axis=1).reset_index(drop=True)
        # Drop the features if the user wants to
        if drop:
            for feature in feature_list:
                self.dataframe.drop(feature, 1)

    # Method to get training and testing datasets
    def _split_dataset(self, scaling=False, test_size=0.2, get_calibration_set=False):
        if get_calibration_set:
            self._extract_calibration_set()
        self.y = self.dataframe["SOLUTION"]
        self.X = self.dataframe.drop("SOLUTION", 1)
        self.X = self.X.drop("PHSR_LOCAL_CC", 1)
        self.X = self.X.drop("PHSR_OVERALL_CC", 1)
        self.X = self.X.drop("RFMC_LOCAL_CC", 1)
        self.X = self.X.drop("RFMC_OVERALL_CC", 1)
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, shuffle=True,
                                                                                test_size=test_size, random_state=42,
                                                                                stratify=self.y)
        self.X_train_unscaled = self.X_train
        self.X_test_unscaled = self.X_test
        self.is_split = True
        if scaling:
            self._scale()

    # Method to balance the number of samples in each class to be 50/50
    def _balance_dataset(self, accross_targets=True):
        self.is_balanced = True
        # Do this accross targets
        if accross_targets:
            new_df = []
            for target in set([x for x in self.dataframe.TARGET_ID]):
                target_data = self.dataframe[self.dataframe.TARGET_ID == target]
                target_solved = target_data[target_data.SOLUTION == 1]
                target_non_solved = target_data[target_data.SOLUTION == 0]
                target_non_solved = target_non_solved.sample(frac=1).reset_index(drop=True)
                new_df.append(
                    pd.concat([target_solved, target_non_solved[:target_solved.shape[0]]]).reset_index(drop=True))
            self.dataframe = pd.concat(new_df)

        # Otherwise ignore the target
        else:
            solved = self.dataframe[self.dataframe.SOLUTION == 1]
            non_solved = self.dataframe[self.dataframe.SOLUTION != 1]
            non_solved = non_solved.sample(frac=1).reset_index(drop=True)
            self.dataframe = pd.concat([solved, non_solved[:solved.shape[0]]]).reset_index(drop=True)

    # Function to scale X
    def _scale(self):
        # Scaling for SVC
        scaler = StandardScaler()
        scaler.fit(self.X_train)
        self.X_train_unscaled = self.X_train
        self.X_test_unscaled = self.X_test
        self.X_train = scaler.transform(self.X_train)
        self.X_test = scaler.transform(self.X_test)
        if self.X_calibration is not None:
            self.X_calibration_unscaled = self.X_calibration
            self.X_calibration = scaler.transform(self.X_calibration)
        self.is_scaled = True

    # Function to get the log of a feature
    def _log_transform(self, colname, new_colname="LOG", insitu=False, force_positive=False):
        # Make a new field
        if not insitu:
            self.dataframe[new_colname] = self.dataframe[colname]
            colname = new_colname
        if force_positive:
            self.dataframe[colname] = self.dataframe[colname].abs()
        # Change values
        self.dataframe[colname] = np.log10(self.dataframe[colname])

    # Function to extract a calibration set
    def _extract_calibration_set(self, fraction_size=0.05):
        if not self.is_balanced:
            self._balance_dataset()
        self.dataframe = self.dataframe.reset_index(drop=True)
        calibration_set = self.dataframe.sample(frac=fraction_size, random_state=22)
        self.dataframe = self.dataframe.drop(calibration_set.index)
        self.y_calibration = calibration_set["SOLUTION"]
        self.X_calibration = calibration_set.drop("SOLUTION", 1)
        self.X_calibration = self.X_calibration.drop("PHSR_LOCAL_CC", 1)
        self.X_calibration = self.X_calibration.drop("PHSR_OVERALL_CC", 1)
        self.X_calibration = self.X_calibration.drop("RFMC_LOCAL_CC", 1)
        self.X_calibration = self.X_calibration.drop("RFMC_OVERALL_CC", 1)
