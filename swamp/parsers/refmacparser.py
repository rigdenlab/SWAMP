from swamp.parsers.parser import Parser


class RefmacParser(Parser):
    """Refmac output parser

    :example:

    >>> from swamp.parsers import RefmacParser
    >>> my_parser = RefmacParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self._rfactor = "NA"
        self._rfree = "NA"
        self._rfactor_delta = ("NA", "NA")
        self._rfree_delta = ("NA", "NA")
        self._bondlenght_delta = ("NA", "NA")
        self._bondangle_delta = ("NA", "NA")
        self._chirvol_delta = ("NA", "NA")

        super(RefmacParser, self).__init__(fname=None, stdout=stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.rfactor, self.rfree, self.rfactor_delta, self.rfree_delta, \
               self.bondlenght_delta, self.bondangle_delta, self.chirvol_delta

    @property
    def rfactor(self):
        return self._rfactor

    @rfactor.setter
    def rfactor(self, value):
        self._rfactor = value

    @property
    def rfactor_delta(self):
        return self._rfactor_delta

    @rfactor_delta.setter
    def rfactor_delta(self, value):
        self._rfactor_delta = value

    @property
    def rfree(self):
        return self._rfree

    @rfree.setter
    def rfree(self, value):
        self._rfree = value

    @property
    def rfree_delta(self):
        return self._rfree_delta

    @rfree_delta.setter
    def rfree_delta(self, value):
        self._rfree_delta = value

    @property
    def bondlenght_delta(self):
        return self._bondlenght_delta

    @bondlenght_delta.setter
    def bondlenght_delta(self, value):
        self._bondlenght_delta = value

    @property
    def bondangle_delta(self):
        return self._bondangle_delta

    @bondangle_delta.setter
    def bondangle_delta(self, value):
        self._bondangle_delta = value

    @property
    def chirvol_delta(self):
        return self._chirvol_delta

    @chirvol_delta.setter
    def chirvol_delta(self, value):
        self._chirvol_delta = value

    def parse(self):
        """Extract the figures of merit from the stdout"""

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
