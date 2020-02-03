import os
import math
from Bio.PDB import *
from swamp.searchmodel_prepare.prepare import PrepareSearchModel
from swamp.searchmodel_prepare.core import Core


class Bfactor(PrepareSearchModel):
    """Wrapper to refactor the bfactors in the search model

    Examples
    --------
    work in progress...

    """

    # ------------------ Class specific properties ------------------

    @property
    def cmd(self):
        """Property to store the command to be executed"""
        return None

    @property
    def modification(self):
        """Property to store the modification to be applied"""
        return "bfactor"

    @property
    def _core_workdir(self):
        """ Workdir where the mr is prepared"""
        return os.path.join(self.workdir, "mr")

    @property
    def _tmp_coreout(self):
        """ Temporary file for the output of the mr preparation"""
        return os.path.join(self._tmp_coreout, "core_out.pdb")

    # ------------------ Class specific methods ------------------

    def prepare(self):

        """Method to prepare the search model using polyala"""

        self.make_workdir()
        os.chdir(self.workdir)
        # Get the mr
        my_core = Core(workdir=self._core_workdir, pdbin=self.pdbin, pdbout=self._tmp_coreout)
        my_core.prepare()
        if my_core.error:
            self.logger.error("Something went wrong preparing the ensemble mr...")
            return
        # Refactor the bfactor
        self.refactor(self._tmp_coreout, self.pdbout)
        self.check_output()

    @staticmethod
    def refactor(pdbin, pdbout):

        """"Method to do the actual bfactor refactoring"""

        input_ensemble = PDBParser().get_structure('bfactor_refactor', pdbin)
        bfactor_dict = Bfactor._get_bfactor_dict(input_ensemble, make_gradient=False)
        for model in input_ensemble.get_models():
            for idx, residue in enumerate(model.get_residues()):
                for atom in residue.get_atoms():
                    atom.set_bfactor(bfactor_dict[model.id][idx])

        # Output the structure
        io = PDBIO()
        io.set_structure(input_ensemble)
        io.save(pdbout)

    @staticmethod
    def _get_atom_dict(input_ensemble):

        """Method to get a dictionary with the atoms present in a Bio.PDBParser"""

        atom_dict = {}
        for model in input_ensemble.get_models():
            for idx, residue in enumerate(model.get_residues()):
                for atom in residue.get_atoms():
                    if atom.get_name() == "CA":
                        if idx in atom_dict.keys():
                            atom_dict[idx].append(atom)
                        else:
                            atom_dict[idx] = [atom]
                        break
        return atom_dict

    @staticmethod
    def _get_model_bfactorgradient(model, lower=10, upper=90):

        """Method to return a list with a gradient of bfactors for a given Bio.PDBParser.model"""

        model_lenght = 0
        for residue in model.get_residues():
            model_lenght += 1
        model_lenght = math.ceil(model_lenght / 2)
        bfactor_values = [lower + x * (upper - lower) / model_lenght for x in range(model_lenght)]

        return bfactor_values[::-1] + bfactor_values

    @staticmethod
    def _get_bfactor_dict(input_ensemble, make_gradient=True, independent_models=False):

        """ Get the bfactor  dictionary"""

        bfactor_dict = {}

        # If the user wants a gradient
        if make_gradient:
            for model in input_ensemble.get_models():
                bfactor_dict[model.id] = Bfactor._get_model_bfactorgradient(model)

        # Otherwise we determine the new bfactor based on ensemble deviation
        else:
            # Get the atoms at each position for all models in a dictionary
            atom_dict = Bfactor._get_atom_dict(input_ensemble)
            # Determine the average distance for each residue position in the ensemble
            all_CA = [x for x in atom_dict.keys()]
            all_CA.sort()
            for residue_position in all_CA:
                # Modify the bfactor per residue in each model
                if not independent_models:
                    distance = 0
                    count = 0
                    for idx, residue_A in enumerate(atom_dict[residue_position]):
                        for residue_B in atom_dict[residue_position][idx + 1:]:
                            distance += residue_A - residue_B
                            count += 1
                    avg_distance = distance / count
                    if avg_distance > 0.8:
                        bfact = 80
                    elif avg_distance <= 0.3:
                        bfact = 10
                    else:
                        bfact = (avg_distance - 0.1) * 100
                    print(residue_position, avg_distance, bfact)
                    # Add the corresponding bfactor to the range in all models
                    for model in input_ensemble.get_models():
                        if model.id in bfactor_dict.keys():
                            bfactor_dict[model.id].append(bfact)
                        else:
                            bfactor_dict[model.id] = [bfact]
                # Modify the bfactors of all models simultaneously
                else:
                    for model_idx, model in enumerate(input_ensemble.get_models()):
                        for residue_idx, residue in enumerate(model.get_residues()):
                            for atom in residue.get_atoms():
                                if atom.get_name() == "CA":
                                    distance=0
                                    count=0
                                    for other_model_residue in atom_dict[residue_idx]:
                                        if atom_dict[residue_idx].index(other_model_residue) != residue_idx:
                                            distance += atom - other_model_residue
                                            count += 1
                                    avg_distance = distance / count
                            if avg_distance > 0.8:
                                bfact = 80
                            elif avg_distance <= 0.3:
                                bfact = 10
                            else:
                                bfact = (avg_distance - 0.1) * 100
                            print(residue_position, avg_distance, bfact)
                            if model.id in bfactor_dict.keys():
                                bfactor_dict[model.id].append(bfact)
                            else:
                                bfactor_dict[model.id] = [bfact]

        return bfactor_dict
