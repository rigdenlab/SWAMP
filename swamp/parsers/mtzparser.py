import re
import gemmi
from enum import Enum
from swamp.parsers.parser import Parser


class MTZLabels(Enum):
    free = re.compile(r"[Ff][Rr][Ee][Ee]")
    i = re.compile(r"I")
    sigi = re.compile(r"SIGI")
    f = re.compile(r"F[P]?$")
    sigf = re.compile(r"SIGF[P]?$")
    i_plus = re.compile(r"I\(\+\)")
    sigi_plus = re.compile(r"SIGI\(\+\)")
    f_plus = re.compile(r"F[P]?\(\+\)")
    sigf_plus = re.compile(r"SIGF[P]?\(\+\)")
    i_minus = re.compile(r"SIGI\(-\)")
    sigi_minus = re.compile(r"SIGI\(-\)")
    f_minus = re.compile(r"F[P]?\(-\)")
    sigf_minus = re.compile(r"SIGF[P]?\(-\)")


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
        """Abstract property to store a summary of the parsed figures of merit"""
        return None

    def parse(self):
        """Method to parse the input mtz file and retrieve tha label names"""

        if self.error:
            self.logger.warning("Previous errors prevent parsing mtz file!")
            return

        for label in MTZLabels:
            matches = list(filter(label.value.match, self.all_labels))
            if any(matches):
                self.__setattr__(label.name, matches[0])
