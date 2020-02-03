import os
import logging
import shutil
from Bio import SeqIO
from pyjob import cexec
from Bio.Alphabet import generic_protein
from swamp.searchmodel_prepare.prepare import PrepareSearchModel


class Molrep(PrepareSearchModel):
    """Molrep wrapper to prepare search model

    Examples
    --------
    work in progress...

    """

    # ------------------ Class specific properties ------------------

    @property
    def cmd(self):
        """Property to store the command to be executed"""
        return "molrep -m {} -s %s -po {}" % self.target_fa

    @property
    def modification(self):
        """Property to store the modification to be applied"""
        return "molrep"

    @property
    def _tmp_fasta(self):
        """Temporary fasta file to be created if necessary"""
        return os.path.join(self.work_dir, "molrep_input_tmp.fasta")

    @property
    def _target_chains(self):
        """Chains present the target fasta file"""
        return list(SeqIO.parse(self.target_fa, "fasta", alphabet=generic_protein))

    @property
    def _target_nchains(self):
        """Number of chains in the target fasta file"""
        return len(self._target_chains)

    @property
    def _molrep_root_template(self):
        """Propety to conain the molrep root template"""
        return os.path.join(self.work_dir, "molrep_out_{}_")

    @property
    def _molrep_align_template(self):
        """Propety to conain the molrep alignment output template"""
        return "%salign.pdb" % self._molrep_root_template

    # ------------------ Class specific methods ------------------

    def make_logfile(self, mode="w"):
        """Method to create the logfile of the wrapper"""

        with open(self.logfile, mode) as fhandle:
            fhandle.write(self.logcontents)

    def _check_fasta_sanity(self):
        """Method to make sure that the target fasta file contains only one chain"""

        if self._target_nchains > 1:
            self.target_fa = self._tmp_fasta
            SeqIO.write(self._target_chains[0], self._tmp_fasta, "fasta")

    def prepare(self):
        """Method to prepare the search model using molrep"""

        self.make_workdir()
        os.chdir(self.workdir)
        self.check_input()
        if self.error:
            return
        self._check_fasta_sanity()

        for model in self.model_list:
            modelID = os.path.basename(model)[:-4]
            cmd = self.cmd.format(model, self._molrep_root_template.format(modelID)).split()
            self.logger.debug(" ".join(cmd))
            self.logcontents = cexec(cmd, stdin=self.keywords, permit_nonzero=True)
            self.make_logfile(mode="a")
            try:
                shutil.copyfile(self._molrep_align_template.format(modelID),
                                self._modified_model_template.format(modelID))
                self.modified_model_list.append(self._modified_model_template.format(modelID))
            except IOError:
                logging.error("Modified model not found! %s" % self._molrep_align_template.format(modelID))
                self.error = True
                pass

        self._merge_models()
        self.check_output()

