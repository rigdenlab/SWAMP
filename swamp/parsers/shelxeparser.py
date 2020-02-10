from swamp.parsers.parser import Parser


class ShelxeParser(Parser):
    """Shelxe output parser

    :example:

    >>> from swamp.parsers.shelxeparser import ShelxeParser
    >>> my_parser = ShelxeParser('<fname>', '<logcontents>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, stdout, logger=None):

        self._cc_eobs_ecalc = "NA"
        self._cc = "NA"
        self._acl = "NA"
        self._average_cc_delta = "NA"
        self._solution = "NO"

        super(ShelxeParser, self).__init__(fname, stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.cc, self.acl, self.cc_eobs_ecalc, self.average_cc_delta, self.solution

    @property
    def cc_eobs_ecalc(self):
        return self._cc_eobs_ecalc

    @cc_eobs_ecalc.setter
    def cc_eobs_ecalc(self, value):
        self._cc_eobs_ecalc = value

    @property
    def cc(self):
        return self._cc

    @cc.setter
    def cc(self, value):
        self._cc = value

    @property
    def acl(self):
        return self._acl

    @acl.setter
    def acl(self, value):
        self._acl = value

    @property
    def solution(self):
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def average_cc_delta(self):
        return self._average_cc_delta

    @average_cc_delta.setter
    def average_cc_delta(self, value):
        self._average_cc_delta = value

    def parse(self):
        """Extract the figures of merit from the logfile and the pdb output file"""

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
