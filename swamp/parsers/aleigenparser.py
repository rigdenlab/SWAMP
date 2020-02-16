from swamp.parsers.parser import Parser


class AleigenParser(Parser):
    """Aleigen contact map alignment output parser

    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar str con_sco: the contact score obtained
    :ivar str c1: number of contacts in the first input contact map
    :ivar str c2: number of contacts in the second input contact map
    :ivar str cmo: contact maximum overlap between the two input contacr maps
    :ivar dict alignemnt: a dictionary with the residue sequence number equivalence across the two input maps

    :example:

    >>> from swamp.parsers import AleigenParser
    >>> my_parser = AleigenParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self.alignment = {}
        self.con_sco = None
        self.c1 = None
        self.c2 = None
        self.cmo = None

        super(AleigenParser, self).__init__(stdout=stdout, fname=None, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.alignment, self.con_sco, self.cmo, self.c1, self.c2

    def parse(self):
        """Extract the figures of merit out of :py:attr:`~swamp.parsers.parser.stdout`"""

        for line in self.stdout.split("\n"):
            line = line.split()
            if len(line) == 4 and line[0] != "Score":
                self.con_sco = float(line[0])
                self.c1 = int(line[1])
                self.c2 = int(line[2])
                self.cmo = float(line[3])
            elif len(line) == 2:
                # residue_b => residue_a
                self.alignment[int(line[1])] = int(line[0])
