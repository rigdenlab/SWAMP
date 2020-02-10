from swamp.parsers.parser import Parser


class AleigenParser(Parser):
    """Aleigen contact map alignment output parser

    :example:

    >>> from swamp.parsers.aleigenparser import AleigenParser
    >>> my_parser = AleigenParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self._alignment = {}
        self._con_sco = None
        self._c1 = None
        self._c2 = None
        self._cmo = None

        super(AleigenParser, self).__init__(stdout=stdout, fname=None, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.alignment, self.con_sco, self.cmo, self.c1, self.c2

    @property
    def c2(self):
        """Property to store c2"""
        return self._c2

    @c2.setter
    def c2(self, value):
        self._c2 = value

    @property
    def cmo(self):
        """Property to store cmo"""
        return self._cmo

    @cmo.setter
    def cmo(self, value):
        self._cmo = value

    @property
    def c1(self):
        """Property to store c1"""
        return self._c1

    @c1.setter
    def c1(self, value):
        self._c1 = value

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

    def parse(self):
        """Extract the figures of merit from the stdout produced by aleigen"""

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
