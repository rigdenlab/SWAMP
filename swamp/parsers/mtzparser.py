import re
import gemmi
from enum import Enum
from swamp.parsers.parser import Parser


class MTZLabels(Enum):
    """An enumerator that contains the regular expression used to detect the column labels of a given MTZ file"""

    free = re.compile(r"[Ff][Rr][Ee][Ee]")
    i = re.compile(r"[Ii]$")
    sigi = re.compile(r"[Ss][Ii][Gg][Ii]$")
    f = re.compile(r"[Ff][Pp]?$")
    sigf = re.compile(r"[Ss][Ii][Gg][Ff][Pp]?$")
    i_plus = re.compile(r"[Ii]\(\+\)")
    sigi_plus = re.compile(r"[Ss][Ii][Gg][Ii]\(\+\)")
    f_plus = re.compile(r"[Ff][Pp]?\(\+\)")
    sigf_plus = re.compile(r"[Ss][Ii][Gg][Ff][Pp]?\(\+\)")
    i_minus = re.compile(r"[Ii]\(-\)")
    sigi_minus = re.compile(r"[Ss][Ii][Gg][Ii]\(-\)")
    f_minus = re.compile(r"[Ff][Pp]?\(-\)")
    sigf_minus = re.compile(r"[Ss][Ii][Gg][Ff][Pp]?\(-\)")


class MtzParser(Parser):
    """Class to parse and store mtz label data.

    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)

    :example:

    >>> from swamp.parsers import MtzParser
    >>> my_parser = MtzParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, logger=None):

        self.reflection_file = None
        self.all_labels = None
        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.free = None
        self.f_plus = None
        self.sigf_plus = None
        self.i_plus = None
        self.sigi_plus = None
        self.f_minus = None
        self.sigf_minus = None
        self.i_minus = None
        self.sigi_minus = None
        self.read_reflections()

        super(MtzParser, self).__init__(fname, logger=logger)

    @property
    def summary(self):
        """Tuple with all the parsed label names"""
        return (self.f, self.sigf, self.i, self.sigi, self.free, self.f_plus, self.sigf_plus, self.i_plus,
                self.sigi_plus, self.f_minus, self.sigf_minus, self.i_minus, self.sigi_minus)

    def read_reflections(self):
        """Read the data in :py:attr:`~swamp.parsers.mtzparser.Mtzparser.reflection_file` file using \
        :py:func:`gemmi.read_mtz_file`"""

        self.reflection_file = gemmi.read_mtz_file(self.fname)
        self.all_labels = [column.label for column in self.reflection_file.columns]

    def parse(self):
        """Parse the input mtz file and retrieve the column names of the labels as described
           at :py:obj:`~swamp.parsers.mtzparser.MTZLabels`"""

        if self.error:
            self.logger.warning("Previous errors prevent parsing mtz file!")
            return

        for label in MTZLabels:
            matches = list(filter(label.value.match, self.all_labels))
            if any(matches):
                self.__setattr__(label.name, matches[0].encode('utf-8'))

        if not any([label for label in self.summary if label is not None]):
            self.logger.error('Cannot find any column names at %s' % self.fname)
            self.error = True
