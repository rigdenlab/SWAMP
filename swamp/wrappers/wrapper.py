import abc
import os
import shutil
import tempfile
import logging

ABC = abc.ABCMeta('ABC', (object,), {})


class Wrapper(ABC):
    """Abstract class for wrappers

    Implements data structures and methods commonly used in wrappers implemented in SWAMP

    :param workdir: working directory
    :type workdir: str, None
    :param silent_start: if True, the logger will not display the start banner (default False)
    :type silent_start: bool
    :param logger: logging interface to be used (default None)
    :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
    :ivar error: if True an error has occurred along the process
    :type error: bool
    ivar logcontents: stores the contents of the log file or stdout
    :type logcontents: str
    :ivar logfile: file name where the logcontents should be written
    type logfile: str
    """

    def __init__(self, workdir=None, logger=None, silent_start=False):

        self._error = False
        self._logcontents = None
        self._silent_start = silent_start
        self.__filthy_files = []
        self._workdir = workdir

        if logger is None:
            self._logger = logging.getLogger(__name__)
        else:
            self._logger = logger

        if workdir is not None:
            self._logfile = os.path.join(self.workdir, '%s.log' % self.wrapper_name)
        else:
            self._logfile = None

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def run(self):
        """Abstract method to run the wrapper"""
        pass

    @abc.abstractmethod
    def get_scores(self, logfile=None):
        """Abstract method to get the scores for the wrapper"""
        pass

    @property
    @abc.abstractmethod
    def keywords(self):
        """Abstract property to store the keywords to use in the wrapper"""
        pass

    @property
    @abc.abstractmethod
    def cmd(self):
        """Abstract property to store command to run in the terminal (if any)"""
        pass

    @property
    @abc.abstractmethod
    def wrapper_name(self):
        """Abstract property to store the name of the wrapper"""
        pass

    @property
    @abc.abstractmethod
    def summary_results(self):
        """Abstract property to store a summary with the results"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def wrapper_header(self):
        """Wrapper header to be displayed in the logger on start"""

        return """\n**********************************************************************
**********************           %s           ********************
**********************************************************************

""" % self.wrapper_name.upper()

    @property
    def logcontents(self):
        return self._logcontents

    @logcontents.setter
    def logcontents(self, value):
        self._logcontents = value

    @property
    def silent_start(self):
        return self._silent_start

    @silent_start.setter
    def silent_start(self, value):
        self._silent_start = value

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def logfile(self):
        return os.path.join(self.workdir, '%s_out.log' % self.wrapper_name)

    @property
    def pdbout(self):
        """Output pdb filename"""
        return os.path.join(self.workdir, '%s_out.pdb' % self.wrapper_name)

    @property
    def mtzout(self):
        """Output mtz filename"""
        return os.path.join(self.workdir, '%s_out.mtz' % self.wrapper_name)

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def _filthy_files(self):
        """Set of files to be deleted after execution"""
        return self.__filthy_files

    @_filthy_files.setter
    def _filthy_files(self, value):
        self.__filthy_files = value

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    # ------------------ Some general methods ------------------

    def make_workdir(self):
        """Create the workdir for the wrapper"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    def _remove_workdir(self):
        """Method to remove the workdir after running wrapper"""

        if os.path.isdir(self.workdir):
            os.chdir(os.path.dirname(os.path.dirname(self.workdir)))
            shutil.rmtree(self.workdir)

    def make_logfile(self, fname=None):
        """Write the log contents into the logfile

        :param fname: specify a file name to write the log (default None)
        :type fname: None, str
        :returns nothing
        :rtype None
        """

        if fname is not None:
            with open(fname, "w") as fhandle:
                fhandle.write(self.logcontents)
        elif self.logfile is not None:
            with open(self.logfile, "w") as fhandle:
                fhandle.write(self.logcontents)
        else:
            self.logger.warning("Impossible to create a logfile, no file name was set!")

    def _cleanup_files(self):
        """Method to clean up temporary files after the wrapper is executed"""

        for fname in self._filthy_files:
            if fname is not None and os.path.isfile(fname):
                os.remove(fname)

    @staticmethod
    def get_cc(workdir, mtzfile, pdbfile, logger=None):
        """Method to get correlation coefficient between mtz file with phase information and a placed search model

        :param workdir: working directory where phenux.get_cc_mtz_pdb will be executed
        :type workdir: str
        :param mtzfile: mtz file name with the phase information
        :type mtzfile: str
        :param pdbfile: pdb file name with the placed search model
        :type pdbfile: str
        :param logger: logging interface to be used on the calulation
        :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
        :returns tuple containing the local CC (str) and the overall CC (str)
        :rtype tuple
        """

        from swamp.wrappers import PhenixCC
        phenix = PhenixCC(pdbin=pdbfile, workdir=workdir, mtzin=mtzfile, logger=logger)
        phenix.run()
        return phenix.local_CC, phenix.overall_CC

    @staticmethod
    def get_tempfile():
        """Method to get a temporary file name

        :returns temporary file name
        :rtype str
        """

        temp_name = next(tempfile._get_candidate_names())
        return os.path.join(os.environ['CCP4_SCR'], '%s' % temp_name)
