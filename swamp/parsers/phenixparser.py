from swamp.parsers.parser import Parser


class PhenixParser(Parser):
    """A phenix.get_cc stdout parser

    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar overall_CC: the overall correlation coefficient parsed from :py:attr:`~swamp.parsers.parser.stdout`
    :ivar local_CC: the local correlation coefficient parsed from :py:attr:`~swamp.parsers.parser.stdout`

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

    def parse(self):
        """Extract the figures of merit from :py:attr:`~swamp.parsers.parser.stdout`"""

        for line in self.stdout.split('\n'):
            if "overall CC" in line:
                self.overall_CC = line.split(":")[1].rstrip().lstrip()
            if "local CC" in line:
                self.local_CC = line.split(":")[1].rstrip().lstrip()

        if self.overall_CC == "NA" or self.local_CC == "NA":
            self.logger.error("Overall / Local CC cannot be found in phenix.get_cc stdout!")
            self.error = True
