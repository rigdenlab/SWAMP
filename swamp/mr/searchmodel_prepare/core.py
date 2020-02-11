import os
import sys
from swamp.mr.searchmodel_prepare.prepare import PrepareSearchModel
sys.path.extend([os.path.join(os.environ['CCP4'], 'lib', 'py2')])
from mrbump.seq_align.MRBUMP_gesamt import Gesamt as MRBUMP_Gesamt


class Core(PrepareSearchModel):
    """Class to trim an ensemble to the core structural alignment

    :param workdir: working directory where the search model will be prepared
    :type workdir: str
    :param pdbin: input pdb file name
    :type pdbin: str
    :param pdbout: output pdb file name
    :type pdbout: str
    :param logger: logging interface for the instance (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`

    :examples

    >>> from swamp.mr.searchmodel_prepare.core import Core
    >>> my_core = Core('<workdir>', '<pdbin>', '<pdbout>')
    >>> my_core.prepare()

    """

    # ------------------ Class specific properties ------------------

    @property
    def modification(self):
        """Property to store the modification to be applied"""
        return "core"

    @property
    def _tmp_csvfile(self):
        """Temporary csv file with the structural alignment between models of the ensemble"""
        return os.path.join(self.workdir, "csvfile.csv")

    @property
    def _tmp_csvCOREfile(self):
        """Temporary csv core file with the core structural alignment between models of the ensemble"""
        return os.path.join(self.workdir, "csvCOREfile.csv")

    @property
    def _tmp_logfile(self):
        """Temporary log file name"""
        return os.path.join(self.workdir, "core_logfile.log")

    @property
    def _tmp_alnfile(self):
        """Temporary alignment file"""
        return os.path.join(self.workdir, "alnfile.ali")

    @property
    def _tmp_scriptfile(self):
        """Temporary script file to execute gesamt"""
        return os.path.join(self.workdir, "scriptfile.sh")

    @property
    def pdbDict(self):
        """A dictionary with the PDB files"""

        pdbDict = dict([])
        for model in self.model_list:
            if len(os.path.basename(model)[:-4]) >= 6:
                pdbDict[model] = "%s" % os.path.basename(model)[:-4][-1]
            else:
                pdbDict[model] = "*"
        return pdbDict

    # ------------------ Class specific methods ------------------

    def _extract_core(self):
        """Method to extract the core structural alignment of a given ensemble"""

        self.logger.debug('Compute structural alignment using gesamt')
        ensemble_truncator = MRBUMP_Gesamt()
        ensemble_truncator.runGesamt(self.model_list, self.pdbDict, seqin="fake.fasta", outputPDB=self._tmp_pdb,
                                     logfile=self._tmp_logfile, alnfile=self._tmp_alnfile, csvfile=self._tmp_csvfile,
                                     script=self._tmp_scriptfile, debug=False)

        self.logger.debug('Truncate ensemble into core')
        ensemble_truncator.makeGesTruncEnsemble(self._tmp_pdb, self.pdbout, variancePercent=10000, sidechain_level=100,
                                                csvFile=self._tmp_csvfile, truncation_level=1000.0,
                                                csvCOREfile=self._tmp_csvCOREfile)

    def prepare(self):
        """Method to prepare the search model ensemble and retrieve the core alignment"""

        self._make_workdir()
        os.chdir(self.workdir)
        self._extract_core()
        self._check_output()