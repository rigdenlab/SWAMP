import gemmi
from Bio import SeqIO
from swamp.parsers import MtzParser
from Bio.Alphabet import generic_protein
from Bio.SeqUtils import molecular_weight


class TargetData(object):
    """Class to store relevant information about the target structure to be solved through MR

    :param str fasta_fname: target's fasta filename
    :param str mtz_fname: target's mtz filename
    :param str phased_mtz_fname: target's mtz filename containing phases (default: None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    """

    def __init__(self, fasta_fname, mtz_fname, phased_mtz_fname=None, logger=None):

        self.fasta_fname = fasta_fname
        self.mtz_fname = mtz_fname
        self.phased_mtz_fname = phased_mtz_fname
        self.mw = None
        self.use_f = False
        self.resolution = None
        self.nreflections = None
        self.spacegroup_symbol = None
        self.solvent = None
        self.seq_length = None
        self.nreflections = None
        self.spacegroup_symbol = None
        self.spacegroup_name = None
        self.ncopies = None
        self.f = None
        self.sigf = None
        self.i = None
        self.sigi = None
        self.free = None
        self.f_plus = None
        self.sigf_plus = None
        self.i_plus = None
        self.sigi_plus = None
        self.f_minus = None
        self.sigf_minus = None
        self.i_minus = None
        self.sigi_minus = None
        self.logger = logger

        self.get_info()

    def get_info(self):
        """Get all the information required to perform MR on the given target and store it into corresponding attributes
        of this :py:obj:`~swamp.mr.targetdata.TargetData` instance"""

        target_chains = [str(chain.seq) for chain in
                         list(SeqIO.parse(self.fasta_fname, "fasta", alphabet=generic_protein))]
        target_chains = list(set(target_chains))
        self.mw = 0.0
        for seq in target_chains:
            seq = seq.replace("X", "A")
            self.mw += round(molecular_weight(seq, "protein"), 2)
        self.seq_length = sum([len(seq) for seq in target_chains])

        mtz_parser = MtzParser(self.mtz_fname, logger=self.logger)
        self.nreflections = mtz_parser.nreflections
        self.spacegroup_symbol = mtz_parser.spacegroup_symbol
        self.spacegroup_name = mtz_parser.spacegroup_symbol.replace(' ', '')
        self.resolution = mtz_parser.resolution

        mtz_parser.parse()
        if mtz_parser.i is None and mtz_parser.f is not None:
            self.use_f = True
        self.estimate_crystal_content(mtz_parser.reflection_file)
        self.f, self.sigf, self.i, self.sigi, self.free, self.f_plus, self.sigf_plus, self.i_plus, self.sigi_plus, \
        self.f_minus, self.sigf_minus, self.i_minus, self.sigi_minus = mtz_parser.summary

    def estimate_crystal_content(self, reflection_file):
        """Estimate the number of copies and the solvent content of the crystal

        :param `gemmi.Mtz` reflection_file: file with the reflections used to estimate the crystal content
        """

        for ncopies in [1, 2, 3, 4, 5]:

            matthews = reflection_file.cell.volume_per_image() / (self.mw * ncopies)
            protein_fraction = 1. / (6.02214e23 * 1e-24 * 1.35 * matthews)
            self.solvent = round((1 - protein_fraction), 1)

            if round(matthews, 3) <= 3.5:
                self.ncopies = ncopies
                break
