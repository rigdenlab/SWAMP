from swamp.detect_solution.data import Data
from lib.parsers.map_align import MapAlignParser
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pandas.plotting import scatter_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold


# Create a class to contain the exploratory analysis of the data
class ExploratoryAnalysis(object):
    # Constructor method
    def __init__(self, csv_file):
        self.data = Data(csv_file)

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

    # Method to get training and testing datasets
    def split_dataset(self, scaling=False, test_size=0.2):
        self.data._split_dataset(scaling, test_size)

    # Make a fucntion to encode the deltas
    def encode_deltas(self, **kwargs):
        self.data._encode_delta(**kwargs)

    # Make the dataset classes balanced
    def balance_dataset(self):
        self.data._balance_dataset()

    # Make functio to encode spacegroup info
    def encode_spacegroup_dataset(self, **kwargs):
        self.data._encode_spacegroup(**kwargs)

    # Function to scale X
    def scale(self):
        self.data._scale()

    # PCA
    def plot_pca(self, fname=None):
        X_std = StandardScaler().fit_transform(self.data.X)
        pca = PCA(.95)
        principalComponents = pca.fit_transform(X_std)
        principalDf = pd.DataFrame(data=principalComponents, columns=['principal component %s' % x for x in
                                                                      range(1, len(pca.explained_variance_ratio_) + 1)])
        finalDf = pd.concat([principalDf, self.data.y], axis=1)
        # Plot PCA
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('Principal Component 1 (%s)' % round(pca.explained_variance_ratio_[0], 2), fontsize=15)
        ax.set_ylabel('Principal Component 2 (%s)' % round(pca.explained_variance_ratio_[1], 2), fontsize=15)
        ax.set_title('2 component PCA', fontsize=20)
        targets = [1.0, 0.0]
        colors = ['g', 'r']
        for target, color in zip(targets, colors):
            indicesToKeep = finalDf['SOLUTION'] == target
            ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                       finalDf.loc[indicesToKeep, 'principal component 2'], c=color, s=5)
        ax.legend(targets)
        ax.grid()
        if fname is not None:
            fig.savefig(fname)
        else:
            fig.show()

        print(pca.explained_variance_ratio_)

    # Scatter matrix all vs all
    def plot_scatter_matrix(self):
        features = [self.data.X_train.columns]
        for feature in features:
            if self.data.X_train.isnull().any().any():
                X_train = self.data.X_train.dropna(axis=1)
            else:
                X_train = self.data.X_train
        attr = list(X_train)
        scatter_matrix(X_train[attr], figsize=(20, 20))
        plt.show()

    def plot_hist(self):
        '''Plot a histogram for each feature.'''
        for name in self.data.X.columns:
            plt.hist(self.data.X[name], bins=20)
            plt.xlabel(name)
            plt.ylabel('number of counts')
            plt.show()

    def plot_density(self, rank_field="CON_RANK", per_target=False, **kwargs):
        solved = self.data.dataframe[self.data.dataframe.SOLUTION == 1]
        # solved = self.data.dataframe[self.data.dataframe.PHSR_LOCAL_CC > 0.5]
        if not per_target:
            ax = solved[rank_field].plot(kind='density', **kwargs)
            ax.set_xlabel(rank_field)
            ax.set_ylabel("Density of solutions")
            # plt.legend(loc='best')
            plt.show()
        else:
            fig, ax = plt.subplots()
            ax.set_xlabel(rank_field)
            ax.set_ylabel("Density of solutions")
            for key, grp in solved.groupby(['TARGET_ID']):
                ax = grp[rank_field].plot(ax=ax, kind='density', label=key, **kwargs)

            plt.legend(loc='best')
            plt.show()

    # Function to merge with map_align data
    def _load_mapalign_data(self, hits_file):
        map_align_parser = MapAlignParser(hits_file)
        map_align_parser.parse_data()
        map_align_parser.dataframe.columns = ['CLST_ID', 'best_params', 'Score', 'cont_sco', 'gap_sco',
                                              'max_scoA', 'max_scoB', 'tot_scoA', 'tot_scoB', 'E_ali', 'E_thread',
                                              'Neff', 'Nali', 'lenA', 'lenB', 'lenAM', 'custom_sco']
        map_align_parser._replace("centroid", "cluster")
        merged_df = self.data.dataframe.merge(map_align_parser.dataframe, how="inner", on="CLST_ID")
        merged_df["Score"] = merged_df["Score"].astype(float)
        merged_df["cont_sco"] = merged_df["cont_sco"].astype(float)
        merged_df["gap_sco"] = merged_df["gap_sco"].astype(float)
        merged_df["custom_sco"] = merged_df["custom_sco"].astype(float)
        self.data.dataframe = merged_df
        del merged_df
        del map_align_parser

    # Function to filter the table so that only the best placement for each model is present
    def get_best_placement(self):
        representative_data = []
        # Get the clst data
        for clst_id in set(self.data.dataframe["CLST_ID"]):
            clst_data = self.data.dataframe[self.data.dataframe.CLST_ID == clst_id]
            # Get the best placement
            clst_data = clst_data.sort_values(by="PHSR_OVERALL_CC", ascending=False)
            clst_representative = clst_data.iloc[[0]]
            representative_data.append(clst_representative)
        self.data.dataframe = pd.concat(representative_data)


if __name__ == '__main__':
    my_analysis = ExploratoryAnalysis(
        "/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_prepared_data_NEWGEN_unbalanced.csv")
    print(my_analysis.data._info_shape)
    my_analysis.omit_samples(feature="MODIF", omit_value="singlechain")

    my_analysis.encode_deltas(initial_field="INIT_RFACT", final_field="FINAL_RFACT", new_fieldname="RFACT_DELTA")
    my_analysis.encode_deltas(initial_field="INIT_RFREE", final_field="FINAL_RFREE", new_fieldname="RFREE_DELTA")
    my_analysis.encode_deltas(initial_field="INIT_CHIRVOL", final_field="FINAL_CHIRVOL", new_fieldname="CHIRVOL_DELTA")
    my_analysis.encode_deltas(initial_field="INIT_BONDLENGHT", final_field="FINAL_BONDLENGHT",
                              new_fieldname="BONDLENGHT_DELTA")
    my_analysis.encode_deltas(initial_field="INIT_BONDANGLE", final_field="FINAL_BONDANGLE",
                              new_fieldname="BONDLANGLE_DELTA")
    my_analysis.data._encode_spacegroup(insitu=True)
    my_analysis.data._log_transform("SHXE_ACL", "log_SHXE_ACL")
    my_analysis.data._log_transform("eLLG", "log_eLLG")
    my_analysis.data._log_transform("TFZ", "log_TFZ")
    my_analysis.data._log_transform("RFZ", "log_RFZ")
    #my_analysis.omit_features(["INIT_RFREE", "INIT_BONDLENGHT", "INIT_BONDANGLE", "INIT_CHIRVOL"])

    my_analysis.omit_features(
        ["TARGET_ID", "CLST_ID", "CON_SCO", "MODIF", 'RESOLUTION', 'RESIDUES_ASU', 'SPACEGROUP', "CHIRVOL_DELTA", "INIT_CHIRVOL", "RFACT_DELTA",
         "FINAL_CHIRVOL", "INIT_RFREE", "FINAL_RFREE", "RFREE_DELTA", 'SHXE_CC', 'SHXE_ACL', "INIT_BONDLENGHT", "FINAL_RFACT",
         "FINAL_BONDLENGHT", "BONDLENGHT_DELTA", "INIT_BONDANGLE", "FINAL_BONDANGLE", "BONDLANGLE_DELTA", "INIT_RFACT",
         'SHXE_Eobs_Ecalc', 'AVG_CC_DELTA', 'N_MODELS', 'CLST_QSCORE', 'CON_RANK', "TOTAL_NRES", "N_RELATED_CLST"])
    # my_analysis.plot_density(bw_method=0.5, xticks=[x for x in range(0, 800, 100)], xlim=(0, 800),
    #                         title="Solutions distribution for all targets")
    # my_analysis.plot_density(per_target=True, bw_method=0.5, xticks=[x for x in range(0, 900, 100)], xlim=(0, 900),
    #                         title="Per target")
    my_analysis.split_dataset()
    my_analysis.plot_scatter_matrix()
    print(my_analysis.data.X.columns)
    my_analysis.plot_pca()
