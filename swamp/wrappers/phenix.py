import os
import shutil
from pyjob import cexec
from swamp.wrappers.wrapper import Wrapper


class PhenixCC(Wrapper):
    """Wrapper around phenix.get_cc_mtz_pdb

     :param workdir: working directory
     :type workdir: str, None
     :param mtzin: file name of the mtz with the phase information to calculate correlation coefficient
     :type mtzin: str
     :param pdbin: pdb file name to calculate correlation coefficient
     :type pdbin: str
     :param logger: logging interface to be used (default None)
     :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
     :param silent_start: if True, the logger will not display the start banner (default False)
     :type silent_start: bool
     :ivar error: if True an error has occurred along the process
     :type error: bool

    :example

    >>> from swamp.wrappers.phenix import PhenixCC
    >>> my_pehnix = PhenixCC('<workdir>', '<pdbin>', '<mtzin>')
    >>> my_pehnix.run()
    """

    def __init__(self, workdir, pdbin, mtzin, logger=None, silent_start=False):

        super(PhenixCC, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self._overall_CC = "NA"
        self._local_CC = "NA"
        self._pdbin = pdbin
        self._mtzin = mtzin

    @property
    def wrapper_name(self):
        return "phenix_cc"

    @property
    def cmd(self):
        """Command to be executed on the shell"""
        return ["phenix.get_cc_mtz_pdb", self.mtzin, self.pdbin]

    @property
    def keywords(self):
        """No keywords are used through stdin in phenix.get_cc_mtz_pdb"""
        return None

    @property
    def summary_results(self):
        return "Results: Local CC - %s   Overall CC - %s" % (self.local_CC, self.overall_CC)

    @property
    def mtzin(self):
        return self._mtzin

    @mtzin.setter
    def mtzin(self, value):
        self._mtzin = value

    @property
    def pdbin(self):
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def local_CC(self):
        """Local correlation coefficient between the mtz file and the pdb file"""
        return self._local_CC

    @local_CC.setter
    def local_CC(self, value):
        self._local_CC = value

    @property
    def overall_CC(self):
        """Overall correlation coefficient between the mtz file and the pdb file"""
        return self._overall_CC

    @overall_CC.setter
    def overall_CC(self, value):
        self._overall_CC = value

    def _check_error(self):
        """Check if errors have occurred during execution"""

        if not os.path.isfile(self.logfile):
            self.logger.error("No phenix output, aborting!")
            self.error = True

    def _clean_files(self):
        """Remove annoying files"""

        try:
            shutil.rmtree(os.path.join(self.workdir, "temp_dir"))
            os.remove(os.path.join(self.workdir, "cc.log"))
            os.remove(os.path.join(self.workdir, "offset.pdb"))
        except OSError:
            self.logger.warning("Unable to clean files after phenix.get_cc! Current wd : %s" % self.workdir)

    def run(self):
        """Run the phenix.get_cc_mtz_pdb"""

        # Manage workdirs
        self.make_workdir()
        current_workdir = os.getcwd()
        os.chdir(self.workdir)

        # Run phenix.get_cc
        self.logger.debug(" ".join(self.cmd))
        self.logcontents = cexec(self.cmd, permit_nonzero=True)
        self.make_logfile()
        self._clean_files()
        self._check_error()
        if self.error:
            return
        self.get_scores()

        # Return to original working directory
        os.chdir(current_workdir)

    def get_scores(self, logfile=None):
        """Parse the log file and obtain the local and the overall CC values

        :param logfile: None
        :returns nothing
        :rtype None
        """

        with open(self.logfile, "r") as fhandle:
            for line in fhandle:
                if "overall CC" in line:
                    self.overall_CC = line.split(":")[1].rstrip().lstrip()
                if "local CC" in line:
                    self.local_CC = line.split(":")[1].rstrip().lstrip()

        if self.overall_CC == "NA" or self.local_CC == "NA":
            self.logger.error("Overall / Local CC not found in phenixCC logfile")
            self.error = True
