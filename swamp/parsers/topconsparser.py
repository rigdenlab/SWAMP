import pandas as pd
from swamp.parsers.parser import Parser
from itertools import groupby
from operator import itemgetter


class TopconsParser(Parser):
    """Topcons file parser

    :param str fname: the file name to be parsed (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the parser (default None)
    :ivar `~pandas.Dataframe` tmhelices: Dataframe with the mappings for the TM helices predicted by topcons
    :ivar `~pandas.Dataframe` residue_topology: Dataframe with the residue topology as predicted by topcons

    :example:

    >>> from swamp.parsers import TopconsParser
    >>> my_parser = TopconsParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname, logger=None):
        self.tmhelices = None
        self.residue_topology = None
        super(TopconsParser, self).__init__(fname, logger=logger)

    @property
    def summary(self):
        """Abstract property to store a summary of the parsed figures of merit"""
        return None

    def parse(self):
        """Parse the :py:attr:`~swamp.parsers.parser.fname` prediction file and retrieve the TM topology"""

        if self.error:
            self.logger.warning("Previous errors prevent parsing TOPCONS file!")
            return

        with open(self.fname, "r") as fhandle:
            self.inputfile_contents = fhandle.readlines()
            try:
                topcons_prediction = self.inputfile_contents[
                    self.inputfile_contents.index('TOPCONS predicted topology:\n') + 1].rstrip()
            except ValueError as e:
                raise ValueError('TOPCONS file cannot be parsed. Please check it is TOPCONS format!')
            residues = []
            for index, ss_residue in enumerate(topcons_prediction):
                if ss_residue == "o":
                    residues.append([index + 1, True, False, False])
                elif ss_residue == "M":
                    residues.append([index + 1, False, True, False])
                else:
                    residues.append([index + 1, False, False, True])

        self.residue_topology = pd.DataFrame(residues)
        self.residue_topology.columns = ["idx", "out", "membr", "in"]
        if any(self.residue_topology.membr.tolist()):
            self._get_tmhelices_map()
        else:
            self.logger.warning('No TM helices were parsed from TOPCONS file!')

    def _get_tmhelices_map(self):
        """Create a datframe with the start/stop of the TM helices contained in the input file"""

        helices = []
        id = 1
        for k, g in groupby(enumerate([x for x in self.residue_topology[self.residue_topology.membr].idx]),
                            lambda ix: ix[0] - ix[1]):
            residues = [x for x in map(itemgetter(1), g)]
            start = residues[0]
            try:
                stop = residues[-1]
            except IndexError:
                stop = start
            helices.append([id, start, stop])
            id += 1
        self.tmhelices = pd.DataFrame(helices)
        self.tmhelices.columns = ["id", "start", "stop"]
