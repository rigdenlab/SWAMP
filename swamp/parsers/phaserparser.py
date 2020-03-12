from swamp.parsers.parser import Parser
from enum import Enum


class PhaserScores(Enum):
    """An enumerator that contains the figures of merit to be parsed from phaser logfile"""

    LLG = 'LLG'
    TFZ = 'TFZ'
    RFZ = 'RFZ'


class PhaserParser(Parser):
    """A Phaser output parser

    :param str fname: the file name to be parsed (default None)
    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar str RFZ: the RFZ as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar str TFZ: the RFZ as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar str LLG: the LLG as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar str eLLG: the eLLG as parsed from the :py:attr:`~swamp.parsers.parser.stdout`
    :ivar str VRMS: the VRMS as parsed from the :py:attr:`~swamp.parsers.parser.stdout`

    :example:

    >>> from swamp.parsers import PhaserParser
    >>> my_parser = PhaserParser('<fname>', '<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, stdout, logger=None):

        self.RFZ = "NA"
        self.TFZ = "NA"
        self.LLG = "NA"
        self.eLLG = "NA"
        self.VRMS = "NA"

        super(PhaserParser, self).__init__(fname, stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the parsed figures of  merit"""
        return self.LLG, self.TFZ, self.RFZ, self.eLLG, self.VRMS

    def parse(self):
        """Extract the figures of merit from :py:attr:`~swamp.parsers.parser.fname` and \
        :py:attr:`~swamp.parsers.parser.stdout`"""

        # Parse the pdbout for LLG, TFZ and RFZ
        with open(self.fname, "r") as fhandle:
            lines = fhandle.readlines()

        llg_remark = [x for x in lines if 'REMARK' in x and "Log-Likelihood Gain" in x]
        figures_of_merit_remark = [x for x in lines if ('REMARK' in x) and ("TFZ" in x or "RFZ" in x)]

        if not any(figures_of_merit_remark):
            self.error = True
            self.logger.error('Cannot find REMARK entry with figures of merit!')
            return

        for attribute in PhaserScores:
            values = [x for x in figures_of_merit_remark[0].split() if '%s=' % attribute.value in x]
            if any(values):
                self.__setattr__(attribute.value, values[-1].split("=")[-1].rstrip().lstrip())

        if any(llg_remark):
            self.LLG = llg_remark[0].split()[-1].rstrip().lstrip()

        # Parse the logfile for eLLG and VRMS
        ellg_reached = False
        for line in self.stdout.split("\n"):
            line = line.rstrip().lstrip()
            if "eLLG   RMSD frac-scat  Ensemble" in line:
                ellg_reached = True
            elif ellg_reached:
                self.eLLG = line.split()[0]
                ellg_reached = False
            if "SOLU ENSEMBLE" in line and "VRMS DELTA" in line:
                self.VRMS = line.split()[5].rstrip().lstrip()
                break

        if self.LLG == "NA" or self.TFZ == "NA":
            self.logger.error("Unable to parse TFZ (%s) and LLG (%s)" % (self.TFZ, self.LLG))
            self.error = True
