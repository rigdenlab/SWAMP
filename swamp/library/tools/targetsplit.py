import os
import gemmi
import itertools
import conkit.io
from swamp.parsers import TopconsParser
from swamp.logger import SwampLogger
from swamp.library.tools.pdb_tools import extract_hierarchy_seqnumber
from swamp.library.tools.contact_tools import extract_interhelical_cmap


class TargetSplit(object):
    """Class to split a given target protein into pairs of TM helices and their interhelical contact information

    :param str workdir: working directory
    :param str conpred: contact prediction file name for the target
    :param str sspred: target's secondary structure prediction file name
    :param str conformat: format of the contact prediction file (default: 'psicov')
    :param str pdb_benchmark: target's pdb file to be used for benchmarking purposes (default: None)
    :ivar subtargets: list of subtargets that the target has been splitted into
    :type subtargets: list

    :example

    >>> from swamp.library.tools.targetsplit import TargetSplit
    >>> splitter = TargetSplit('<workdir>', '<conpred>', '<sspred>')
    >>> splitter.split()

    """

    def __init__(self, workdir, conpred, sspred, conformat="psicov", pdb_benchmark=None, logger=None):
        self._workdir = workdir
        self._conpred = conkit.io.read(conpred, conformat).top_map
        self._make_workdir()
        self._subtargets = None
        self._subtargets_pdb = None
        self._pdb_benchmark = pdb_benchmark
        self._error = False
        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(logfile=None, use_console=True, console_level='info')
        else:
            self._logger = logger
        self._sspred = TopconsParser(sspred, logger=self.logger)
        self.sspred.parse()
        if self.sspred.error:
            self.logger.warning('Previous errors detected while parsing TM topology prediction!')
            self.error = True

    # ------------------ Some general properties ------------------

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def pdb_benchmark(self):
        return self._pdb_benchmark

    @pdb_benchmark.setter
    def pdb_benchmark(self, value):
        self._pdb_benchmark = value

    @property
    def subtargets(self):
        return self._subtargets

    @subtargets.setter
    def subtargets(self, value):
        self._subtargets = value

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def conformatsspred(self):
        return self._conformatsspred

    @conformatsspred.setter
    def conformatsspred(self, value):
        self._conformatsspred = value

    @property
    def conpred(self):
        return self._conpred

    @conpred.setter
    def conpred(self, value):
        self._conpred = value

    @property
    def sspred(self):
        return self._sspred

    @sspred.setter
    def sspred(self, value):
        self._sspred = value

    @property
    def _possible_helical_pairs(self):
        """List of all possible helical pair combinations"""
        return list(itertools.combinations(self.sspred.tmhelices.id, 2))

    @property
    def ranked_subtargets(self):
        """A list of the subtargets sorted according to the number of interhelical contacts"""

        if self.subtargets is None:
            return None
        else:
            ranked_scores = [cmap.ncontacts for cmap in self.subtargets]
            ranked_scores.sort()
            ranked_scores.reverse()
            ranked_cmaps = [-1 for x in ranked_scores]
            # Loop through the targets
            for cmap in self.subtargets:
                ranked_cmaps[ranked_scores.index(cmap.ncontacts)] = cmap
                ranked_scores[ranked_scores.index(cmap.ncontacts)] = -1
            # Eliminate targets with score 0
            ranked_cmaps = [cmap for cmap in ranked_cmaps if cmap.ncontacts > 0]
            # Store the ranked targets
            return ranked_cmaps

    @property
    def subtargets_pdb(self):
        """List containing the :obj:`gemmi.Structure` of the subtargets (only if pdb_benchmark is not None)"""

        if self.pdb_benchmark is None or self.subtargets is None:
            return None
        elif self._subtargets_pdb is None:
            self.extract_pdb_tmhelices(self.pdb_benchmark)
        return self._subtargets_pdb

    # ------------------ Some general methods ------------------

    def _make_workdir(self):
        """Make the working directory"""

        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def _helix_range(self, helix_id):
        """Extract the range of residue numbers that a given tm helix spans.

        :param helix_id: the id of the tm helix of interest
        :type helix_id: str, int
        :returns a list with the range of residues that the tm helix of interest spans
        :rtype list
        """

        if self.sspred is not None and self.sspred.tmhelices is not None:
            return range(int(self.sspred.tmhelices[self.sspred.tmhelices.id == helix_id].start),
                         int(self.sspred.tmhelices[self.sspred.tmhelices.id == helix_id].stop) + 1)
        else:
            return None

    def get_helical_pair(self, helical_pair):
        """Method to extract the residues spanning a given helical pair

        :param helical_pair: indicates the pair of helices to extract
        :type: helical_pair: tuple, list
        :returns region_A, region_B: residues sequence numbers of the regions where the tm helices of the pair span
        :rtype tuple
        """

        region_A = [x for x in self._helix_range(helical_pair[0])]
        region_B = [x for x in self._helix_range(helical_pair[1])]

        return region_A, region_B

    def get_interhelical_contacts(self, helical_pair, score_threshold=0.45):
        """Get the interhelical contacts beyond certain score threshold of a helical pair

        :param helical_pair: a tuple indicatinh the helical pair of interest
        :type helical_pair: tuple
        :param score_threshold: the score threshold to include a contact in the resulting contact map
        :type score_threshold: float
        :returns a contact map instance with the interhelical contacts that met the threshold
        :rtype `:obj:conkit.core.ContactMap`
        """

        helix_region_A, helix_region_B = self.get_helical_pair(helical_pair)
        helices = (helix_region_A, helix_region_B)
        residues = helix_region_A + helix_region_B
        residues.sort()
        new_id = "_".join(map(str, helical_pair))

        return extract_interhelical_cmap(self.conpred, helices, residues, new_id, score_threshold=score_threshold)

    def extract_pdb_tmhelices(self, pdbfile):
        """Extract the regions of the :obj:`gemmi.Structure` that correspond with the each of the tm helices
         predicted in the secondary structure file

         :param pdbfile: location of the pdbfile of the target structure
         :type pdbfile: str
         :returns nothing
         :rtype None
         """

        self._subtargets_pdb = {}
        hierarchy = gemmi.read_structure(pdbfile)
        hierarchy.remove_ligands_and_waters()
        for sbtrgt in self.subtargets:
            tmhelices = list(map(int, sbtrgt.id.split("_")))
            region = [inner for outer in self.get_helical_pair(tmhelices) for inner in outer]
            self._subtargets_pdb[sbtrgt.id] = extract_hierarchy_seqnumber(hierarchy, seq_numbers=region)

    def split(self):
        """Split the target into its contacting helical pairs. It will probe all possible combinations
        of helical pairs."""

        if self.sspred.error:
            self.logger.warning('Previous errors prevented splitting the target into helical pairs!')
            self.error = True
            return

        self.subtargets = []
        for helical_pair in self._possible_helical_pairs:
            interhelical_cmap = self.get_interhelical_contacts(helical_pair)
            self.subtargets.append(interhelical_cmap)
        self.logger.info("Splitted the target into %s helical pairs" % len(self.subtargets))
