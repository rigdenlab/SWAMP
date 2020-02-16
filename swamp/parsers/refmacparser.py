from swamp.parsers.parser import Parser


class RefmacParser(Parser):
    """Refmac output parser

    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar str rfactor: the Rfactor as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar tuple rfactor_delta: the intial and final Rfactor as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar str rfree: the Rfree as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar tuple rfree_delta: the intial and final Rfree as parsed from the :py:attr:`~swamp.parsers.parser.fname`
    :ivar tuple bondlenght_delta: the intial and final Bond length parsed from :py:attr:`~swamp.parsers.parser.fname`
    :ivar tuple bondangle_delta: the intial and final Bond angle parsed from :py:attr:`~swamp.parsers.parser.fname`
    :ivar tuple chirvol_delta: the intial and final Chirvol as parsed from the :py:attr:`~swamp.parsers.parser.fname`

    :example:

    >>> from swamp.parsers import RefmacParser
    >>> my_parser = RefmacParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self.rfactor = "NA"
        self.rfree = "NA"
        self.rfactor_delta = ("NA", "NA")
        self.rfree_delta = ("NA", "NA")
        self.bondlenght_delta = ("NA", "NA")
        self.bondangle_delta = ("NA", "NA")
        self.chirvol_delta = ("NA", "NA")

        super(RefmacParser, self).__init__(fname=None, stdout=stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.rfactor, self.rfree, self.rfactor_delta, self.rfree_delta, \
               self.bondlenght_delta, self.bondangle_delta, self.chirvol_delta

    def parse(self):
        """Extract the figures of merit from :py:attr:`~swamp.parsers.parser.stdout`"""

        reached_end = False
        for line in self.stdout.split("\n"):
            if "Final results" in line:
                reached_end = True
            elif reached_end and "R factor" in line:
                self.rfactor = line.split()[3].rstrip().encode('utf-8')
                self.rfactor_delta = (line.split()[2].rstrip().encode('utf-8'), self.rfactor)
            elif reached_end and "R free" in line:
                self.rfree = line.split()[3].rstrip().encode('utf-8')
                self.rfree_delta = (line.split()[2].rstrip().encode('utf-8'), self.rfree)
            elif reached_end and "Rms BondLength" in line:
                self.bondlenght_delta = (
                    line.split()[2].rstrip().encode('utf-8'), line.split()[3].rstrip().encode('utf-8'))
            elif reached_end and "Rms BondAngle" in line:
                self.bondangle_delta = (
                    line.split()[2].rstrip().encode('utf-8'), line.split()[3].rstrip().encode('utf-8'))
            elif reached_end and "Rms ChirVolume" in line:
                self.chirvol_delta = (
                    line.split()[2].rstrip().encode('utf-8'), line.split()[3].rstrip().encode('utf-8'))

        # If there is no rfree or rfactor, there was an error
        if self.rfactor == "NA" and self.rfree == "NA":
            self.logger.error("Refmac did not report Rfree and Rfactor !")
            self.error = True
