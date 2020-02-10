import abc
import os
import logging

ABC = abc.ABCMeta('ABC', (object,), {})


class Parser(ABC):
    """Parser abstract class

    This class contains general methods and data structures to extract information from output and prediction files

    :argument fname: name of the file to be parsed (default None)
    :type fname: str, None
    :argument stdout: the stdout to be parsed (default None)
    :type str, None
    :ivar error: True if an error has occurred while parsing the file
    :type error: bool
    :ivar _inputfile_contents: contains the full contents of the input file
    :type _inputfile_contents: str, None
    """

    def __init__(self, fname=None, stdout=None, logger=None):
        self._fname = fname
        self._error = False
        self._stdout = stdout
        self._inputfile_contents = None
        self._check_input()
        if logger is None:
            self._logger = logging.getLogger(__name__)
        else:
            self._logger = logger

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def parse(self):
        """ Abstract method to run the parser"""
        pass

    @abc.abstractmethod
    @property
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def stdout(self):
        return self._stdout

    @stdout.setter
    def stdout(self, value):
        self._stdout = value

    @property
    def inputfile_contents(self):
        return self._inputfile_contents

    @inputfile_contents.setter
    def inputfile_contents(self, value):
        self._inputfile_contents = value

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        self._fname = value

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    # ------------------ Some general methods ------------------

    def _check_input(self):
        """Check if the input file name exists"""

        if self.fname is not None and not os.path.isfile(self.fname):
            self.error = True
            self.logger.error('Cannot find input file, please make sure it exists: %s' % self.fname)
