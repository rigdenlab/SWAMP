import os
import sys
from swamp.searchmodel_prepare.prepare import PrepareSearchModel
sys.path.extend([os.path.join(os.environ['CCP4'], 'lib', 'py2')])
from mrbump.seq_align.MRBUMP_gesamt import Gesamt as MRBUMP_Gesamt


class Core(PrepareSearchModel):
    """Core wrapper to prepare search model

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

        return "polyALA"

    @property
    def _tmp_csvfile(self):

        """Temporary csv file"""

        return os.path.join(self.workdir, "csvfile.csv")

    @property
    def _tmp_csvCOREfile(self):

        """Temporary csv mr file"""

        return os.path.join(self.workdir, "csvCOREfile.csv")

    @property
    def _tmp_logfile(self):

        """Temporary log file"""

        return os.path.join(self.workdir, "core_logfile.log")

    @property
    def _tmp_alnfile(self):

        """Temporary alignment file"""

        return os.path.join(self.workdir, "alnfile.ali")

    @property
    def _tmp_scriptfile(self):

        """Temporary script file"""

        return os.path.join(self.workdir, "scriptfile.sh")

    @property
    def pdbDict(self):

        """Property a dictionary with the PDB files (sth related with molrep)"""

        pdbDict = dict([])
        for model in self.model_list:
            if len(os.path.basename(model)[:-4]) >= 6:
                pdbDict[model] = "%s" % os.path.basename(model)[:-4][-1]
            else:
                pdbDict[model] = "*"
        return pdbDict

    # ------------------ Class specific methods ------------------

    def _extract_core(self):

        """Method to extract the mr of a given ensemble"""

        ensemble_truncator = MRBUMP_Gesamt()
        ensemble_truncator.runGesamt(self.model_list, self.pdbDict, seqin="fake.fasta", outputPDB=self._tmp_pdb,
                                     logfile=self._tmp_logfile, alnfile=self._tmp_alnfile, csvfile=self._tmp_csvfile,
                                     script=self._tmp_scriptfile, debug=False)
        ensemble_truncator.makeGesTruncEnsemble(self._tmp_pdb, self.pdbout, variancePercent=10000, sidechain_level=100,
                                                csvFile=self._tmp_csvfile, truncation_level=1000.0,
                                                csvCOREfile=self._tmp_csvCOREfile)

    def prepare(self):

        """Method to prepare the mr ensemble of the search model"""

        self.make_workdir()
        os.chdir(self.workdir)
        self._extract_core()
        self.check_output()
