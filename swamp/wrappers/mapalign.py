import os
import swamp
import conkit.io
import numpy as np
from pyjob import cexec
from swamp.wrappers.wrapper import Wrapper
from swamp.wrappers.gesamt import Gesamt


class MapAlign(Wrapper):
    """Wrapper around map_align

    :param workdir: working directory
    :type workdir: str
    :param map_a: file name of the map A to be used in the alignment
    :type map_a: str
    :param pdb_a: file name of the pdb A to be used for benchmarking of the alignment
    :type pdb_a: str
    :param format_a: contact prediction file format of map A (default 'mapalign')
    :type format_a: str
    :param map_b: file name of the map B to be used in the alignment
    :type map_b: str
    :param pdb_b: file name of the pdb B to be used for benchmarking of the alignment
    :type pdb_b: str
    :param format_b: contact prediction file format of map B (default 'mapalign')
    :type format_b: str
    :param gap_o: gap opening penalty (default -1)
    :type gap_o: int
    :param gap_e: gap extension penalty (default -0.01)
    :type gap_e: float
    :param sep_cut: separation cutoff to compute contacts (default 0)
    :type sep_cut: int
    :param logger: logging interface for the wrapper (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`
    :ivar error: if True an error has occurred along the process
    :type error: bool

    :examples

    >>> from swamp.wrappers.mapalign import MapAlign
    >>> mapalign = MapAlign('<workdir>', '<map_a>', '<map_b>')
    >>> mapalign.run()

    """

    def __init__(self, workdir, map_a, map_b, format_a='mapalign', format_b='mapalign', gap_o=-1, gap_e=-0.01,
                 sep_cut=0, pdb_a=None, pdb_b=None, logger=None):

        super(MapAlign, self).__init__(workdir=workdir, logger=logger)

        self._alignment = None
        self._con_sco = None
        self._gap_sco = None
        self._total_sco = None
        self._alignment_length = None
        self._sep_cut = sep_cut
        self._gap_e = gap_e
        self._gap_o = gap_o
        self._qscore = np.nan
        self._rmsd = np.nan
        self._seq_id = np.nan
        self._n_align = np.nan
        self._pdb_a = pdb_a
        self._format_a = format_a
        self._format_b = format_b
        self._map_a = map_a
        self._map_b = map_b
        self._pdb_b = pdb_b

    # ------------------ Some properties ------------------

    @property
    def wrapper_name(self):
        return "mapalign"

    @property
    def source(self):
        """Location of the executable file of al-eigen"""
        return swamp.SRC_MAPALIGN

    @property
    def cmd(self):
        """Command to be executed in the shell"""
        return [self.source] + self.keywords

    @property
    def input_a(self):
        """File name of the map file A to use in the alignment"""
        if self.format_a != self._reference_mapformat:
            return os.path.join(self.workdir, "{}.{}").format(self.map_a_id, self._reference_mapformat)
        else:
            return self.map_a

    @property
    def input_b(self):
        """File name of the map file B to use in the alignment"""
        if self.format_b != self._reference_mapformat:
            return os.path.join(self.workdir, "{}.{}").format(self.map_b_id, self._reference_mapformat)
        else:
            return self.map_b

    @property
    def keywords(self):
        """Keywords used to calculate the alignment"""
        return ['-a', self.input_a, '-b', self.input_b, "-sep_cut", str(self.sep_cut), "-gap_e", str(self.gap_e),
                "-gap_o", str(self.gap_o), '-silent']

    @property
    def summary_results(self):
        """List with the figures of merit of the alignment"""
        return [self.map_a, self.map_b, self.con_sco, self.gap_sco, self.total_sco, self.alignment_length, self.qscore,
                self.rmsd, self.seq_id, self.n_align]

    @property
    def map_a(self):
        return self._map_a

    @map_a.setter
    def map_a(self, value):
        self._map_a = value

    @property
    def format_a(self):
        return self._format_a

    @format_a.setter
    def format_a(self, value):
        self._format_a = value

    @property
    def format_b(self):
        return self._format_b

    @format_b.setter
    def format_b(self, value):
        self._format_b = value

    @property
    def map_b(self):
        return self._map_b

    @map_b.setter
    def map_b(self, value):
        self._map_b = value

    @property
    def pdb_a(self):
        return self._pdb_a

    @pdb_a.setter
    def pdb_a(self, value):
        self._pdb_a = value

    @property
    def pdb_b(self):
        return self._pdb_b

    @pdb_b.setter
    def pdb_b(self, value):
        self._pdb_b = value

    @property
    def sep_cut(self):
        return self._sep_cut

    @sep_cut.setter
    def sep_cut(self, value):
        self._sep_cut = value

    @property
    def gap_o(self):
        return self._gap_o

    @gap_o.setter
    def gap_o(self, value):
        self._gap_o = value

    @property
    def gap_e(self):
        return self._gap_e

    @gap_e.setter
    def gap_e(self, value):
        self._gap_e = value

    @property
    def rmsd(self):
        """Property to store rmsd of the structural alignment (only if pdb benchmark is not None)"""
        return self._rmsd

    @rmsd.setter
    def rmsd(self, value):
        self._rmsd = value

    @property
    def n_align(self):
        """Property to store no. of aligned residues of the structural alignment (only if pdb benchmark is not None)"""
        return self._n_align

    @n_align.setter
    def n_align(self, value):
        self._n_align = value

    @property
    def seq_id(self):
        """Property to store the sequence identity of the structural alignment (only if pdb benchmark is not None)"""
        return self._seq_id

    @seq_id.setter
    def seq_id(self, value):
        self._seq_id = value

    @property
    def qscore(self):
        """Property to store qscore of the structural alignment (only if pdb benchmark is not None)"""
        return self._qscore

    @qscore.setter
    def qscore(self, value):
        self._qscore = value

    @property
    def map_a_id(self):
        """Basename of the map A file without extension"""
        return os.path.splitext(os.path.basename(self.map_a))[0]

    @property
    def map_b_id(self):
        """Basename of the map B file without extension"""
        return os.path.splitext(os.path.basename(self.map_b))[0]

    @property
    def alignment(self):
        """Property to store  (dictionary with the residue sequence number equivalence across the input maps)"""
        return self._alignment

    @alignment.setter
    def alignment(self, value):
        self._alignment = value

    @property
    def con_sco(self):
        """Property to store the contact score obtained with the map alignment"""
        return self._con_sco

    @con_sco.setter
    def con_sco(self, value):
        self._con_sco = value

    @property
    def gap_sco(self):
        """Property to store gap score assigned to the map alignment"""
        return self._gap_sco

    @gap_sco.setter
    def gap_sco(self, value):
        self._gap_sco = value

    @property
    def total_sco(self):
        """Property to store total score (contact_score - gap_score)"""
        return self._total_sco

    @total_sco.setter
    def total_sco(self, value):
        self._total_sco = value

    @property
    def alignment_length(self):
        """Property to store contact map alignment length"""
        return self._alignment_length

    @alignment_length.setter
    def alignment_length(self, value):
        self._alignment_length = value

    @property
    def _reference_mapformat(self):
        """The reference contact map format to be used with the aligner"""
        return 'mapalign'

    # ------------------ Some general methods ------------------

    def _check_alignment_input(self):
        """Check that input maps are in correct format, if not convert them"""

        if self.format_a != self._reference_mapformat:
            self.make_workdir()
            self.logger.debug((self.map_a, self.format_a, self.input_a, self._reference_mapformat))
            self._filthy_files.append(self.input_a)
            conkit.io.convert(self.map_a, self.format_a, self.input_a, self._reference_mapformat)
        if self.format_b != self._reference_mapformat:
            self.make_workdir()
            self.logger.debug((self.map_b, self.format_b, self.input_b, self._reference_mapformat))
            self._filthy_files.append(self.input_b)
            conkit.io.convert(self.map_b, self.format_b, self.input_b, self._reference_mapformat)

    def get_scores(self, logfile=None):
        """Method to parse the stdout and extract the scores and residue alignment

        :param logfile: Not in use
        :type logfile: None
        :returns nothing
        :rtype None
        """

        self.alignment = {}
        for line in self.logcontents.split("\n"):
            if line != "" and line.split()[0] == "MAX":
                line = line.split()
                self.con_sco = float(line[4])
                self.gap_sco = float(line[5])
                self.total_sco = float(line[6])
                self.alignment_length = float(line[7])
                for residue in line[8:]:
                    residue = residue.split(":")
                    # residue_b => residue_a
                    self.alignment[int(residue[1])] = int(residue[0])

    def benchmark(self):
        """Perform a structural alignment between the actual structures for benchmarking purposes"""

        gesamt = Gesamt.get_optimum_alignment((self.pdb_a, self.pdb_b))[0]
        if gesamt.error:
            self.rmsd = np.nan
            self.seq_id = np.nan
            self.qscore = np.nan
            self.n_align = np.nan
        else:
            self.rmsd = gesamt.rmsd
            self.seq_id = gesamt.seq_id
            self.qscore = gesamt.qscore
            self.n_align = gesamt.n_align

    def run(self):
        """Obtain the optimal alignment between both maps

        :returns nothing
        :rtype None
        :raises EnvironmentError if the executable binary file of the aligner is not found in the system
        """

        if self.source is None:
            raise EnvironmentError("Cannot find %s binary executable!" % self.wrapper_name)

        self._check_alignment_input()
        self.logger.debug(" ".join(self.cmd))
        self.logcontents = cexec(self.cmd, permit_nonzero=True)
        if self.logcontents == b'':
            self.error = True
            self.logger.error("No output found for the contact map alignment!")
            return
        self.logger.debug(self.logcontents)

        self.get_scores()

        if self.pdb_a is not None and self.pdb_b is not None:
            self.benchmark()

        self._cleanup_files()
