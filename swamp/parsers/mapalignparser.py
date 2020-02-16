from swamp.parsers.parser import Parser


class MapAlignParser(Parser):
    """map_align contact map alignment output parser

    :example:

    >>> from swamp.parsers import MapAlignParser
    >>> my_parser = MapAlignParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self._alignment = {}
        self._con_sco = None
        self._gap_sco = None
        self._total_sco = None
        self._alignment_length = None

        super(MapAlignParser, self).__init__(stdout=stdout, fname=None, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.alignment, self.con_sco, self.gap_sco, self.total_sco, self.alignment_length

    @property
    def alignment(self):
        """Property to store  (dictionary with the residue sequence number equivalence across the input maps)"""
        return self._alignment

    @alignment.setter
    def alignment(self, value):
        self._alignment = value

    @property
    def con_sco(self):
        """Property to store the contact score obtained with the map alignment"""
        return self._con_sco

    @con_sco.setter
    def con_sco(self, value):
        self._con_sco = value

    @property
    def gap_sco(self):
        """Property to store gap score assigned to the map alignment"""
        return self._gap_sco

    @gap_sco.setter
    def gap_sco(self, value):
        self._gap_sco = value

    @property
    def total_sco(self):
        """Property to store total score (contact_score - gap_score)"""
        return self._total_sco

    @total_sco.setter
    def total_sco(self, value):
        self._total_sco = value

    @property
    def alignment_length(self):
        """Property to store contact map alignment length"""
        return self._alignment_length

    @alignment_length.setter
    def alignment_length(self, value):
        self._alignment_length = value

    def parse(self):
        """Extract the figures of merit from the stdout produced by map_align"""

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
