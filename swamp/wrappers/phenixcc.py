import os
import shutil
from pyjob import cexec
from swamp.wrappers.wrapper import Wrapper
from swamp.parsers import PhenixParser


class PhenixCC(Wrapper):
    """Wrapper around phenix.get_cc_mtz_pdb

    :param str workdir: working directory
    :param str mtzin: file name of the mtz with the phase information to calculate correlation coefficient
    :param str pdbin: pdb file name to calculate correlation coefficient
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :param bool silent_start: if True, the logger will not display the start banner (default False)
    :ivar bool error: if True an error has occurred along the process
    :ivar str ovarall_CC: overall correlation coefficient between the phases in the mtz file and the pdb file
    :ivar str local_CC: local correlation coefficient between the phases in the mtz file and the pdb file

    :example:

    >>> from swamp.wrappers import PhenixCC
    >>> my_pehnix = PhenixCC('<workdir>', '<pdbin>', '<mtzin>')
    >>> my_pehnix.run()
    """

    def __init__(self, workdir, pdbin, mtzin, logger=None, silent_start=False):

        super(PhenixCC, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self.overall_CC = "NA"
        self.local_CC = "NA"
        self.pdbin = pdbin
        self.mtzin = mtzin

    @property
    def wrapper_name(self):
        """The name of this `~swamp.wrapper.wrapper.Wrapper` child class (phenix_cc)"""
        return "phenix.get_cc_mtz_pdb"

    @property
    def source(self):
        """Override :py:attr:`~swamp.wrappers.wrapper.Wrapper.source` to include PHENIX/build/bin instead of CCP4/bin"""
        return os.path.join(os.environ['PHENIX'], 'build', 'bin', self.wrapper_name)

    @property
    def cmd(self):
        """Command to be executed on the shell"""
        return [self.source, self.mtzin, self.pdbin]

    @property
    def keywords(self):
        """No keywords are used through stdin in phenix.get_cc_mtz_pdb"""
        return None

    @property
    def summary_results(self):
        """A summary with the obtained figures of merit"""
        return "Results: Local CC - %s   Overall CC - %s" % (self.local_CC, self.overall_CC)

    def _check_error(self):
        """Check if errors have occurred during execution of :py:attr:`~swamp.wrappers.phenixcc.PhenixCC.cmd`"""

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

    def _run(self):
        """Run :py:attr:`~swamp.wrappers.phenixcc.PhenixCC.cmd` and store the results"""

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
        """Use :py:obj:`~swamp.parsers.phenixparser.PhenixParser` to parse the log file and obtain values for \
        :py:attr:`~swamp.wrappers.phenixcc.PhenixCC.local_CC` and \
        :py:attr:`~swamp.wrappers.phenixcc.PhenixCC.overall_CC`

        :param logfile: Not in use
        """

        parser = PhenixParser(logger=self.logger, stdout=self.logcontents)
        parser.parse()

        if parser.error:
            self.error = True
            self.logger.warning('Previous errors while parsing phenix.get_cc output detected!')
        else:
            self.local_CC, self.overall_CC = parser.summary
