from swamp.parsers.parser import Parser


class PhenixParser(Parser):
    """A phenix.get_cc output parser

    :example:

    >>> from swamp.parsers import PhenixParser
    >>> my_parser = PhenixParser('<stdout>')
    >>> my_parser.parse()
    """

    def __init__(self, stdout, logger=None):

        self._overall_CC = "NA"
        self._local_CC = "NA"

        super(PhenixParser, self).__init__(stdout=stdout, fname=None, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.local_CC, self.overall_CC

    @property
    def local_CC(self):
        """Local correlation coefficient between the mtz file and the pdb file"""
        return self._local_CC

    @local_CC.setter
    def local_CC(self, value):
        self._local_CC = value

    @property
    def overall_CC(self):
        """Overall correlation coefficient between the mtz file and the pdb file"""
        return self._overall_CC

    @overall_CC.setter
    def overall_CC(self, value):
        self._overall_CC = value

    def parse(self):
        """Extract the figures of merit from the logfile and the pdb output file"""

        for line in self.stdout.split('\n'):
            if "overall CC" in line:
                self.overall_CC = line.split(":")[1].rstrip().lstrip()
            if "local CC" in line:
                self.local_CC = line.split(":")[1].rstrip().lstrip()

        if self.overall_CC == "NA" or self.local_CC == "NA":
            self.logger.error("Overall / Local CC cannot be found in phenix.get_cc stdout!")
            self.error = True
