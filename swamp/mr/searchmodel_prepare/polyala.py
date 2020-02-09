import os
import gemmi
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from swamp.library.tools.pdb_tools import renumber_hierarchy
from swamp.mr.searchmodel_prepare.prepare import PrepareSearchModel


class PolyALA(PrepareSearchModel):
    """PolyALA wrapper to prepare search model

    Examples
    --------
    >>> from swamp.searchmodel_prepare.polyala import PolyALA
    >>> my_polyala = PolyALA('<workdir>', '<pdbin>', '<pdbout>')
    >>> my_polyala.prepare()

    """

    # ------------------ Class specific properties ------------------

    @property
    def cmd(self):
        """Property to store the command to be executed"""
        return None

    @property
    def modification(self):
        """Property to store the modification to be applied"""
        return "polyALA"

    @property
    def _target_chains(self):
        """Chains present the target fasta file"""
        return list(SeqIO.parse(self.target_fa, "fasta", alphabet=generic_protein))

    @property
    def _target_nchains(self):
        """Number of chains in the target fasta file"""
        return len(self._target_chains)

    # ------------------ Class specific methods ------------------

    def prepare(self):

        """Method to prepare the search model using polyala"""

        self.make_workdir()
        os.chdir(self.workdir)

        for model in self.model_list:
            modelID = os.path.basename(model)[:-4]
            modified_model = self._modified_model_template.format(modelID)
            self.truncate_polyALA(pdbin=model, pdbout=modified_model)
            self.transfer_flags_pdb(pdb_ref=model, pdb_file=modified_model)
            self.modified_model_list.append(modified_model)

        self._merge_models()
        self.check_output()

    @staticmethod
    def truncate_polyALA(pdbin, pdbout):

        """Method to truncate a given pdb into polyala"""

        original_hierarchy = gemmi.read_structure(pdbin)
        original_hierarchy.remove_ligands_and_waters()

        for residue in original_hierarchy[0][0]:
            residue.trim_to_alanine()
            residue.name = "ALA"

        renumber_hierarchy(original_hierarchy)

        original_hierarchy.write_pdb(pdbout)

    @staticmethod
    def transfer_flags_pdb(pdb_ref, pdb_file, flags_to_transfer=("CRYST1", "SCALE", "REMARK"), overwrite=True):

        """ Method to transfer PDB flags between two given pdb files"""

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
