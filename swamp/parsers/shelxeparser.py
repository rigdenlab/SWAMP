from swamp.parsers.parser import Parser


class ShelxeParser(Parser):
    """Shelxe output parser

    :param str fname: the file name to be parsed (default None)
    :param str stdout: the stdout to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar str cc_eobs_ecalc: the correlation coeff. between Eobs and Ecalc
    :ivar str cc: correlation coeff. obtained with the best tracing cycle
    :ivar str acl: average chain length obtained with the best tracing cycle
    :ivar str average_cc_delta: the average delta of correlation coeff. between tracing cycles
    :ivar str solution: 'YES' if correlation coeff > 25 otherwise 'NO'

    :example:

    >>> from swamp.parsers import ShelxeParser
    >>> my_parser = ShelxeParser('<fname>', '<logcontents>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, stdout, logger=None):

        self.cc_eobs_ecalc = "NA"
        self.cc = "NA"
        self.acl = "NA"
        self.average_cc_delta = "NA"
        self.solution = "NO"

        super(ShelxeParser, self).__init__(fname, stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.cc, self.acl, self.cc_eobs_ecalc, self.average_cc_delta, self.solution

    def parse(self):
        """Extract the figures of merit from :py:attr:`~swamp.parsers.parser.fname` and \
        :py:attr:`~swamp.parsers.parser.stdout`"""

        # Get cc average delta from the stdout
        cc_scores = []
        for line in self.stdout.split("\n"):
            if 'Overall CC between native Eobs and Ecalc (from fragment)' in line:
                self.cc_eobs_ecalc = line.split("=")[1].lstrip().rstrip()[:-1]
            elif 'CC for partial structure against native data =' in line:
                cc_scores.append(
                    float(line.rstrip().lstrip().split("=")[1].lstrip().replace(" %", "").lstrip().rstrip()))
        self.average_cc_delta = sum([cc_scores[n] - cc_scores[n - 1] for n in range(1, len(cc_scores))]) / (
                len(cc_scores) - 1)

        # Get the cc and acl from the output file
        with open(self.fname, "r") as fhandle:
            for line in fhandle:
                if "TITLE" in line:
                    self.cc = line.split("=")[1].split("%")[0].rstrip().lstrip()
                    shelxe_residues = float(line.split("%")[1].split()[0].rstrip().lstrip())
                    shelxe_chains = int(line.split("%")[1].split()[3].rstrip().lstrip())
                    self.acl = str(round(shelxe_residues / shelxe_chains))
                    break
        fhandle.close()

        # Determine if there was a solution
        if self.cc != "NA" and float(self.cc) >= 25:
            self.solution = "YES"
