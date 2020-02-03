import pandas as pd
from swamp.parsers.parser import Parser
from itertools import groupby
from operator import itemgetter


class TopconsParser(Parser):
    """Topcons file parser

    :ivar tmhelices: Dataframe with the mappings for the TM helices predicted by topcons
    :type tmhelices: :object pandas.Dataframe, None
    :ivar residue_topology: Dataframe with the residue topology as predicted by topcons
    :type residue_topology: :object pandas.Dataframe, None

    :example:

    >>> from swamp.parsers.topcons import TopconsParser
    >>> my_parser = TopconsParser('<fname>')
    >>> my_parser.parse()
    """

    def __init__(self, fname):
        self._tmhelices = None
        self._residue_topology = None
        super(TopconsParser, self).__init__(fname)

    @property
    def tmhelices(self):
        return self._tmhelices

    @tmhelices.setter
    def tmhelices(self, value):
        self._tmhelices = value

    @property
    def residue_topology(self):
        return self._residue_topology

    @residue_topology.setter
    def residue_topology(self, value):
        self._residue_topology = value

    def parse(self):
        """Parse the topcons prediction file and retrieve the TM topology"""

        self._check_input()
        if self.error:
            self.logger.error("Error detected, please check input is correct...")
            return

        with open(self.fname, "r") as fhandle:
            self.inputfile_contents = fhandle.readlines()
            topcons_prediction = self.inputfile_contents[0].rstrip()
            residues = []
            for index, ss_residue in enumerate(topcons_prediction):
                if ss_residue == "i":
                    residues.append([index + 1, True, False, False])
                elif ss_residue == "M":
                    residues.append([index + 1, False, True, False])
                else:
                    residues.append([index + 1, False, False, True])

            self.residue_topology = pd.DataFrame(residues)
            self.residue_topology.columns = ["idx", "out", "membr", "in"]

        self._get_tmhelices_map()

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
