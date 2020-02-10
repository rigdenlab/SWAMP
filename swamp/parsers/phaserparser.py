from swamp.parsers.parser import Parser


class PhaserParser(Parser):
    """Phaser output parser

    :example:

    >>> from swamp.parsers.phaserparser import PhaserParser
    >>> my_parser = PhaserParser('<fname>', '<logcontents>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, stdout, logger=None):

        self._RFZ = "NA"
        self._TFZ = "NA"
        self._LLG = "NA"
        self._eLLG = "NA"
        self._VRMS = "NA"

        super(PhaserParser, self).__init__(fname, stdout, logger=logger)

    @property
    def summary(self):
        """A summary of the figures of  merit"""
        return self.LLG, self.TFZ, self.RFZ, self.eLLG, self.VRMS

    @property
    def RFZ(self):
        return self._RFZ

    @RFZ.setter
    def RFZ(self, value):
        self._RFZ = value

    @property
    def TFZ(self):
        return self._TFZ

    @TFZ.setter
    def TFZ(self, value):
        self._TFZ = value

    @property
    def LLG(self):
        return self._LLG

    @LLG.setter
    def LLG(self, value):
        self._LLG = value

    @property
    def eLLG(self):
        return self._eLLG

    @eLLG.setter
    def eLLG(self, value):
        self._eLLG = value

    @property
    def VRMS(self):
        return self._VRMS

    @VRMS.setter
    def VRMS(self, value):
        self._VRMS = value

    def parse(self):
        """Extract the figures of merit from the logfile and the pdb output file"""

        # Parse the pdbout for LLG, TFZ and RFZ
        with open(self.fname, "r") as fhandle:
            for line in fhandle:
                if "Log-Likelihood Gain" in line:
                    self.LLG = line.split(":")[1].rstrip().lstrip()
                    continue
                elif "TFZ" in line or "RFZ" in line:
                    for score in line.split():
                        if "TFZ==" in score:
                            self.TFZ = score.split("==")[1].rstrip().lstrip()
                        elif "TFZ=" in score:
                            self.TFZ = score.split("=")[1].rstrip().lstrip()
                        elif "RFZ==" in score:
                            self.RFZ = score.split("==")[1].rstrip().lstrip()
                        elif "RFZ=" in score:
                            self.RFZ = score.split("=")[1].rstrip().lstrip()
                    continue
                # If LLG and TFZ are stored, break
                if self.TFZ != "NA" and self.LLG != "NA" and self.RFZ != "NA":
                    break

        # Parse the logfile for eLLG and VRMS
        ellg_reached = False
        for line in self.stdout.split("\n"):
            line = line.rstrip().lstrip()
            if "eLLG   RMSD frac-scat  Ensemble" in line:
                ellg_reached = True
            elif ellg_reached:
                self.eLLG = line.split()[0]
                ellg_reached = False
            if "SOLU ENSEMBLE" in line and "VRMS DELTA" in line:
                self.VRMS = line.split()[5].rstrip().lstrip()
                break

        if self.LLG == "NA" or self.TFZ == "NA":
            self.logger.error("Unable to parse TFZ (%s) and LLG (%s)" % (self.TFZ, self.LLG))
            self.error = True
