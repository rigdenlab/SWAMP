# Make the imports
import os
import sys
import collections
import subprocess
import dill as pickle
import mmtbx.model
import iotbx.pdb
import shlex
import iotbx.pdb.amino_acid_codes
from simbad.util.mtz_util import crystal_data, GetLabels
from simbad.mr.anomalous_util import AnodeSearch
from simbad.parsers.anode_parser import AnodeParser
import pandas as pd
import tempfile

package_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

sys.path.append(os.path.join(package_path, "swamp", "wrappers"))
from run_phaser import Phaser
from run_refmac import Refmac
from run_shelxe import Shelxe
from run_crank2 import Crank2

sys.path.append(os.path.join(package_path, "swamp", "searchmodel_prepare"))
from prepare_search_model import PrepareSearchModel, resequence_pdb_file

sys.path.append(os.path.join(package_path, "lib"))
from truncate_polyALA import PolyALATruncator

sys.path.append(os.path.join(package_path, "lib", "parsers"))
from make_customDB import renumber_residues

# Create the necessary named tuples
SearchModel = collections.namedtuple("SearchModelInfo",
                                     ["pdbfile", "eRMSD", "n_search", "modification", "total_nresidues", "clst_id",
                                      "con_sco_norm", "n_models", "clst_RMSD", "con_rank"])
Target = collections.namedtuple("TargetInfo",
                                ["pdbfile", "fastafile", "mtzfile", "cmapfile", "phasedmtz", "solvent", "nchains_asu",
                                 "anomalous_signal", "MW", "resolution", "nresidues_chain", "subtarget_con",
                                 "subtarget_fa", "con_format"])
PhaserInfo = collections.namedtuple("PhaserInfo", ["sgalternative", "early_kill", "workdir", "log", "mtzout", "pdbout",
                                                   "root", "timeout", "disable_check", "xyzout", "threads", "use_tncs"])
RefmacInfo = collections.namedtuple("RefmacInfo", ["workdir", "pdbout", "mtzout", "log", "make_hydr", "ncyc",
                                                   "weight_matrix", "ridg_dist_sigm"])
FixedSolution = collections.namedtuple("FixedSolution",
                                       ["sol_file", "pdbfile", "mtzfile", "eRMSD", "ident", "LLG", "TFZ", "local_CC",
                                        "overall_CC", "rfact", "rfree"])
Crank2Info = collections.namedtuple("Crank2Info",
                                    ["dirout", "xyzout", "hklout", "keyin", "logout", "anomalous_scatterer",
                                     "wavelength_type", "xyzin", "mode", "stdout"])


# Create a class to parse the results of a given swamp run
class MrExtractResults(object):
    # Constructor method
    def __init__(self, directory):
        self.mother_directory = directory
        self.pickle_list = []
        self.results_list = [
            ["TARGET_ID", "CLST_ID", "RESOLUTION", "RESIDUES_ASU", "SPACEGROUP", "LLG", "TFZ", "RFZ", "VRMS", "eLLG",
             "PHSR_LOCAL_CC", "PHSR_OVERALL_CC", "INIT_RFACT", "FINAL_RFACT", "INIT_RFREE", "FINAL_RFREE",
             "INIT_BONDLENGHT", "FINAL_BONDLENGHT", "INIT_BONDANGLE", "FINAL_BONDANGLE", "INIT_CHIRVOL",
             "FINAL_CHIRVOL", "RFMC_LOCAL_CC", "RFMC_OVERALL_CC", "SHXE_CC", "SHXE_ACL", "SHXE_Eobs_Ecalc",
             "AVG_CC_DELTA", "CON_SCO", "N_MODELS", "MODIF", "eRMSD", "CLST_QSCORE", "N_RELATED_CLST", "TOTAL_NRES",
             "CON_RANK", "SOLUTION"]]
        for directory in os.listdir(self.mother_directory):
            if not os.path.isdir(os.path.join(self.mother_directory, directory)):
                continue
            if os.path.isfile(os.path.join(self.mother_directory, directory, "results.pckl")):
                self.pickle_list.append(os.path.join(self.mother_directory, directory, "results.pckl"))

    # Function to unpickle a results file
    @staticmethod
    def _unpickle(pickle_file):
        results = []
        with open(pickle_file, "rb") as pickle_handle:
            for x in range(0, 33):
                try:
                    results.append(pickle.load(pickle_handle))
                except EOFError:
                    break
        pickle_handle.close()
        return results

    # Function to assist to obtain con_rank (sometimes in the old pipeline this would be missing)
    @staticmethod
    def _get_con_ranks(rank_file):
        rank_dictionary = {}
        with open(rank_file, "r") as f_handle:
            cluster_template = "cluster_{}"
            for line in f_handle:
                line = line.rstrip().split()
                rank = line[0]
                if rank == "Rank":
                    continue
                cluster_id = cluster_template.format(line[2].split("_")[1])
                rank_dictionary[cluster_id] = rank
        return rank_dictionary

    # Function to assist to obtain the clst_qscore (at the moment the pipeline has a bug preventing record this)
    @staticmethod
    def _get_clst_qscore(search_model_dir):
        qscore_dict = {}
        search_model_list = [os.path.join(search_model_dir, file) for file in os.listdir(search_model_dir)]
        for searchmodel in search_model_list:
            # get hierarchy
            pdb_inp = iotbx.pdb.input(file_name=searchmodel)
            model = mmtbx.model.manager(model_input=pdb_inp)
            pdb_hierarchy = model.get_hierarchy()
            # Loop through data structures to get the number of models and residues
            n_models = 0
            for model in pdb_hierarchy.models():
                n_models += 1
            # Get the qscore
            gesamtEXE = os.path.join(os.environ["CCP4"], "bin", "gesamt")
            for model in range(1, n_models + 1):
                gesamtEXE += " %s -s /%s/" % (searchmodel, model)
            process_args = shlex.split(gesamtEXE, posix=False)
            p = subprocess.Popen(process_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, stderr = p.communicate()
            avg_qscore = "NA"
            for line in stdout.split("\n"):
                if n_models == 2 and "Q-score" in line:
                    avg_qscore = line.rstrip().lstrip().split(":")[1].split()[0].rstrip().lstrip()
                    break
                elif n_models > 2 and "quality Q" in line:
                    avg_qscore = line.rstrip().lstrip().split(":")[1].split()[0].rstrip().lstrip()
                    break
            if avg_qscore == "NA":
                print(search_model_dir)
                print(gesamtEXE)
                print("Something went wrong getting the avg. Qscore of the search model! %s" % searchmodel)
            qscore_dict[os.path.basename(searchmodel)[:-4]] = avg_qscore
        # Return dictionary
        return qscore_dict

    # Function to assist to get related clusters (this is actally not something we record on the pipeline...)
    @staticmethod
    def _get_related_clst(gesamt_centroid_archive, centroid_dir, nthreads):
        related_clst_dict = {}
        centroid_list = [os.path.join(centroid_dir, file) for file in os.listdir(centroid_dir)]
        gesamtEXE = "%s {} -archive %s -nthreads=%s -min1=0 -min2=0 -o {}" % (
            os.path.join(os.environ["CCP4"], "bin", "gesamt"), gesamt_centroid_archive, nthreads)
        for centroid in centroid_list:
            # Get output file
            tmp_file = os.path.join(tempfile._get_default_tempdir(), "%s_hits.txt" % os.path.basename(centroid[:-4]))
            # Scan the archive
            cmd = gesamtEXE.format(centroid, tmp_file)
            process_args = shlex.split(cmd, posix=False)
            p = subprocess.Popen(process_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, stderr = p.communicate()
            related_clst = 0
            with open(tmp_file, "r") as f_handle:
                for line in f_handle:
                    if line[0] != "#":
                        line = line.rstrip().lstrip().split()
                        if float(line[2]) > 0.6 and float(line[2]) != 1.0:
                            related_clst += 1
                        elif float(line[2]) < 0.6:
                            break
            related_clst_dict["cluster_%s" % os.path.basename(centroid).split("_")[1][:-4]] = related_clst
            os.remove(tmp_file)
        # Return dictionary
        return related_clst_dict

    # Function to extract all results
    def extract(self, verbose=False, rank_file=None, search_model_dir=None, centroid_dir=None,
                centroid_gesamt_archive=None, nthreads=1):
        # Get dictionaries with the extra features
        if rank_file is not None:
            rank_dict = self._get_con_ranks(rank_file)
        else:
            rank_dict = None
        if search_model_dir is not None:
            qscore_dict = self._get_clst_qscore(search_model_dir)
        else:
            qscore_dict = None
        if centroid_dir is not None and centroid_gesamt_archive is not None:
            related_dict = self._get_related_clst(gesamt_centroid_archive=centroid_gesamt_archive,
                                                  centroid_dir=centroid_dir, nthreads=nthreads)
        else:
            related_dict = None
        # Extract data from the pickles
        for pickle_file in self.pickle_list:
            if verbose:
                print("\tunpickle: %s" % pickle_file)
            pickled_runs = self._unpickle(pickle_file)
            for run in pickled_runs:
                if verbose:
                    print("\t\trun: %s %s %s" % (
                        run.search_models[0].clst_id, run.search_models[0].modification, run.search_models[0].ermsd))
                current_run_results = []
                current_run_results.append(os.path.basename(run.target.mtzfile)[:4])
                current_run_results.append(run.search_models[0].clst_id)
                current_run_results.append(str(run.target.resolution))
                current_run_results.append(str(int(run.target.nresidues_chain) * int(run.target.nchains_asu)))
                current_run_results.append(run.target.spacegroup_symbol.replace(' ', ''))
                current_run_results.append(run.phaser.LLG)
                current_run_results.append(run.phaser.TFZ)
                current_run_results.append(run.phaser.RFZ)
                current_run_results.append(run.phaser.VRMS)
                current_run_results.append(run.phaser.eLLG)
                current_run_results.append(run.phaser.local_CC)
                current_run_results.append(run.phaser.overall_CC)
                current_run_results.append(run.refmac.rfactor_delta[0])
                current_run_results.append(run.refmac.rfactor_delta[1])
                current_run_results.append(run.refmac.rfree_delta[0])
                current_run_results.append(run.refmac.rfree_delta[1])
                current_run_results.append(run.refmac.bondlenght_delta[0])
                current_run_results.append(run.refmac.bondlenght_delta[1])
                current_run_results.append(run.refmac.bondangle_delta[0])
                current_run_results.append(run.refmac.bondangle_delta[1])
                current_run_results.append(run.refmac.chirvol_delta[0])
                current_run_results.append(run.refmac.chirvol_delta[1])
                current_run_results.append(run.refmac.local_CC)
                current_run_results.append(run.refmac.overall_CC)
                current_run_results.append(run.shelxe.cc)
                current_run_results.append(run.shelxe.acl)
                current_run_results.append(run.shelxe.cc_eobs_ecalc)
                current_run_results.append(run.shelxe.average_cc_delta)
                current_run_results.append(
                    run.search_models[0].con_sco_norm if run.search_models[0].con_sco_norm is not None else "NA")
                current_run_results.append(str(run.search_models[0].n_models))
                current_run_results.append(run.search_models[0].modification)
                current_run_results.append(run.search_models[0].ermsd)
                if qscore_dict is not None:
                    current_run_results.append(qscore_dict[run.search_models[0].clst_id])
                else:
                    current_run_results.append("NA")
                if related_dict is not None:
                    current_run_results.append(related_dict[run.search_models[0].clst_id])
                else:
                    current_run_results.append("NA")
                current_run_results.append(str(run.search_models[0].total_nresidues))
                if rank_dict is not None:
                    current_run_results.append(
                        run.search_models[0].con_rank if run.search_models[0].con_rank is not None else rank_dict[
                            run.search_models[0].clst_id])
                else:
                    current_run_results.append(
                        run.search_models[0].con_rank if run.search_models[0].con_rank is not None else "NA")
                if run.phaser.local_CC != "NA" and float(run.phaser.local_CC) > 0.2:
                    current_run_results.append("1")
                else:
                    current_run_results.append("0")

                self.results_list.append(current_run_results)

    # Method to write results into a file
    def save_file(self, outfile):
        with open(outfile, "w") as f_handle:
            line_template = "{}\n"
            for line in self.results_list:
                f_handle.write(line_template.format("\t".join(line)))

    # Method to return a pandas dataframe with the data
    def get_pd_dataframe(self, make_proportional=False, omit_NA=None):
        my_dframe = pd.DataFrame(self.results_list[1:], columns=self.results_list[0])
        if omit_NA is not None:
            for field in omit_NA:
                my_dframe = my_dframe[my_dframe[str(field)] != "NA"]
        if not make_proportional:
            return my_dframe
        else:
            return self.make_50_50(my_dframe)

    # Method to drop necessary failures to make proportions 50/50
    @staticmethod
    def make_50_50(pandas_df):
        solved = pandas_df[pandas_df.SOLUTION == "1"]
        non_solved = pandas_df[pandas_df.SOLUTION == "0"]
        non_solved = non_solved.sample(frac=1).reset_index(drop=True)
        return pd.concat([solved, non_solved[:solved.shape[0]]]).reset_index(drop=True)


if __name__ == '__main__':
    _OMIT_NA_FIELDS = ["RESOLUTION", "RESIDUES_ASU", "LLG", "TFZ", "RFZ", "VRMS", "eLLG", "FINAL_RFACT", "FINAL_RFREE",
                       "SHXE_CC", "SHXE_ACL", "SHXE_Eobs_Ecalc", "N_MODELS", "TOTAL_NRES", "CON_RANK", "SOLUTION"]

    # Get the saved results
    old_data_filtered = pd.read_csv(
        "/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_filtered_data_unbalanced.csv",
        na_filter=False, skipinitialspace=True, thousands=',')
    old_data_raw = pd.read_csv(
        "/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_filtered_data_unbalanced.csv",
        na_filter=False, skipinitialspace=True, thousands=',')
    new_data_filtered = [old_data_filtered]
    new_data_raw = [old_data_raw]
    # Extract new results
    #new_data_filtered = []
    #new_data_raw = []
    for pdbID in ['4ryo']:
        print("Extracting: %s" % pdbID)
        my_results = MrExtractResults(
            "/mnt/sda1/MR_edge_cases/%s_MR/fragment_picking/MR_results/" % pdbID)
        my_results.extract(verbose=True, nthreads=45,
                           centroid_gesamt_archive="/mnt/sda1/MR_edge_cases/%s_MR/fragment_picking/gesamt_archive" % pdbID,
                           centroid_dir="/mnt/sda1/MR_edge_cases/%s_MR/fragment_picking/centroids" % pdbID,
                           search_model_dir="/mnt/sda1/MR_edge_cases/%s_MR/fragment_picking/ensembles" % pdbID,
                           rank_file="/mnt/sda1/MR_edge_cases/%s_MR/fragment_picking/ranked_hits.list" % pdbID)
        my_data_filtered = my_results.get_pd_dataframe(make_proportional=False, omit_NA=_OMIT_NA_FIELDS)
        new_data_filtered.append(my_data_filtered)
        my_data_raw = my_results.get_pd_dataframe(make_proportional=False)
        new_data_raw.append(my_data_raw)

    # CONCATENATE
    data_filtered = pd.concat(new_data_filtered).reset_index(drop=True)
    data_raw = pd.concat(new_data_raw).reset_index(drop=True)

    # SAVE FILE
    data_filtered.to_csv(
        "/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_filtered_data_unbalanced.csv",
        index=False)
    data_raw.to_csv("/home/filo/PycharmProjects/CON-MOL/swamp/detect_solution/training_data/ml_raw_data_unbalanced.csv",
                    index=False)
