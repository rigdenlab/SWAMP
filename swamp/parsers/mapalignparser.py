from swamp.parsers.parser import Parser


class MapAlignParser(Parser):
    """map_align contact map alignment output parser

    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar str con_sco: the contact score obtained after alignment of the input maps
    :ivar str gap_sco: the gap score obtained in the alignment
    :ivar str total_sco: the gap score obtained after alignment of the input maps
    :ivar int alignment_length: the length of the alignment obtained between the two input maps
    :ivar dict alignemnt: a dictionary with the residue sequence number equivalence across the two input maps

    :example:

    >>> from swamp.parsers import MapAlignParser
    >>> my_parser = MapAlignParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self.alignment = {}
        self.con_sco = None
        self.gap_sco = None
        self.total_sco = None
        self.alignment_length = None

        super(MapAlignParser, self).__init__(stdout=stdout, fname=None, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.alignment, self.con_sco, self.gap_sco, self.total_sco, self.alignment_length

    def parse(self):
        """Extract the figures of merit out of :py:attr:`~swamp.parsers.parser.stdout`"""

        self.alignment = {}
        for line in self.stdout.split("\n"):
            if line != "" and line.split()[0] == "MAX":
                line = line.split()
                self.con_sco = float(line[4])
                self.gap_sco = float(line[5])
                self.total_sco = float(line[6])
                self.alignment_length = float(line[7])
                for residue in line[8:]:
                    residue = residue.split(":")
                    # residue_b => residue_a
                    self.alignment[int(residue[1])] = int(residue[0])
