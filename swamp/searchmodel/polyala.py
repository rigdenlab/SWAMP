import os
import gemmi
from swamp.library.tools.pdb_tools import renumber_hierarchy
from swamp.searchmodel import SearchModel


class PolyALA(SearchModel):
    """Class to strip side chains from the search model.

    Extends :py:obj:`~swamp.searchmodel.searchmodel.SearchModel`

    :examples

    >>> from swamp.searchmodel.polyala import PolyALA
    >>> my_polyala = PolyALA('<workdir>', '<pdbin>', '<pdbout>')
    >>> my_polyala.prepare()

    """

    # ------------------ Class specific properties ------------------

    @property
    def modification(self):
        """Property to store the modification to be applied ("polyala")"""
        return "polyala"

    # ------------------ Class specific methods ------------------

    def prepare(self):
        """Trim the side chains out of the models in the :py:attr:`~swamp.searchmodel.prepare.model_list`"""

        self._make_workdir()
        os.chdir(self.workdir)

        for idx, model in enumerate(self.model_list):
            modelID = os.path.basename(model)[:-4]
            modified_model = self._modified_model_template.format(modelID)
            self.logger.debug('Truncating model %s %s -> %s' % (idx, modelID, modified_model))
            self.truncate_polyALA(pdbin=model, pdbout=modified_model)
            self.logger.debug('Transfer flags to new pdb file')
            self.transfer_flags_pdb(pdb_ref=model, pdb_file=modified_model)
            self.modified_model_list.append(modified_model)

        self.logger.debug('Merge models into ensemble')
        self._merge_models()
        self._check_output()

    @staticmethod
    def truncate_polyALA(pdbin, pdbout):
        """Method to truncate a given pdb into poly-ala

        :param str pdbin: input pdb file name
        :param str pdbout: output pdb file name
        """

        original_hierarchy = gemmi.read_structure(pdbin)
        original_hierarchy.remove_ligands_and_waters()

        for residue in original_hierarchy[0][0]:
            residue.trim_to_alanine()
            residue.name = "ALA"

        renumber_hierarchy(original_hierarchy)

        original_hierarchy.write_pdb(pdbout)

    @staticmethod
    def transfer_flags_pdb(pdb_ref, pdb_file, flags_to_transfer=("CRYST1", "SCALE", "REMARK"), overwrite=True):
        """Transfer PDB flags between two given pdb files

        :param str pdb_ref: pdb file with the reference flags to be transferred
        :param str pdb_file: pdb file where the flags will be transferred
        :param tuple flags_to_transfer: set of flags that need to be transferred
        :param bool overwrite: if False, pdb_file original flags will be  kept (default True)
        """

        with open(pdb_ref, "r") as pdbreference, open(pdb_file, "r") as pdbfile:
            # Read in the lines to transfer
            lines_to_include = []
            for line in pdbreference:
                if line.split()[0] in flags_to_transfer:
                    lines_to_include.append(line)
            # Store the lines of the pdbfile
            for line in pdbfile:
                if not overwrite:
                    lines_to_include.append(line)
                elif line.split()[0] not in flags_to_transfer:
                    lines_to_include.append(line)
        # Re-open the file, this time write mode
        pdbreference.close()
        pdbfile.close()
        with open(pdb_file, "w") as pdbfile:
            for line in lines_to_include:
                pdbfile.write(line)
