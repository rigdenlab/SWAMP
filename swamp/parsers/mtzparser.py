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
    i_minus = re.compile(r"[Ss][Ii][Gg][Ii]\(-\)")
    sigi_minus = re.compile(r"[Ss][Ii][Gg][Ii]\(-\)")
    f_minus = re.compile(r"[Ff][Pp]?\(-\)")
    sigf_minus = re.compile(r"[Ss][Ii][Gg][Ff][Pp]?\(-\)")


class MtzParser(Parser):
    """Class to parse and store mtz label data.

    :example:

    >>> from swamp.parsers.mtzparser import MtzParser
    >>> my_parser = MtzParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, logger=None):

        self.reflection_file = gemmi.read_mtz_file(fname)
        self.all_labels = [column.label for column in self.reflection_file.columns]
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

        super(MtzParser, self).__init__(fname, logger=logger)

    @property
    def summary(self):
        """Abstract property not in use"""
        return None

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
