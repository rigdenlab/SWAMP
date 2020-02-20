import abc
import os
import logging

ABC = abc.ABCMeta('ABC', (object,), {})


class Parser(ABC):
    """Parser abstract class

    This class contains general methods and data structures to extract information from the output created by the \
     classes at :py:obj`~swamp.wrappers`

    :param str fname: name of the file to be parsed (default None)
    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar bool error: True if an error has occurred while parsing the file or the stdout
    :ivar str inputfile_contents: contains the full contents of the input file if any is given
    """

    def __init__(self, fname=None, stdout=None, logger=None):
        self.fname = fname
        self.error = False
        self.stdout = stdout
        self.inputfile_contents = None
        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger
        self.check_input()

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def parse(self):
        """ Abstract method to run the parser"""
        pass

    @property
    @abc.abstractmethod
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit"""
        pass

    # ------------------ Some general methods ------------------

    def check_input(self):
        """Check if :py:attr:`~swamp.parsers.parser.fname` exists"""

        if self.fname is not None and not os.path.isfile(self.fname):
            self.error = True
            self.logger.error('Cannot find input file, please make sure it exists: %s' % self.fname)

