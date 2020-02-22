import os
import conkit.io
import numpy as np
from pyjob import cexec
from swamp.wrappers.wrapper import Wrapper
from swamp.wrappers.gesamt import Gesamt
from swamp.parsers import MapAlignParser


class MapAlign(Wrapper):
    """Wrapper around map_align

    :param str workdir: working directory
    :param str map_a: file name of the map A to be used in the alignment
    :param str pdb_a: file name of the pdb A to be used for benchmarking of the alignment
    :param str format_a: contact prediction file format of map A (default 'mapalign')
    :param str map_b: file name of the map B to be used in the alignment
    :param str pdb_b: file name of the pdb B to be used for benchmarking of the alignment
    :param str format_b: contact prediction file format of map B (default 'mapalign')
    :param int gap_o: gap opening penalty (default -1)
    :param float gap_e: gap extension penalty (default -0.01)
    :param int sep_cut: separation cutoff to compute contacts (default 0)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :ivar bool error: if True an error has occurred along the process
    :ivar str con_sco: the contact score obtained after alignment of the input maps
    :ivar str gap_sco: the gap score obtained in the alignment
    :ivar str total_sco: the gap score obtained after alignment of the input maps
    :ivar int alignment_length: the length of the alignment obtained between the two input maps
    :ivar float qscore: qscore as reported by :py:obj:`~swamp.wrappers.gesamt.Gesamt`
    :ivar float rmsd: the obtained rmsd as reported by :py:obj:`~swamp.wrappers.gesamt.Gesamt`
    :ivar float seq_id: seq identity between the input structures obtained with :py:obj:`~swamp.wrappers.gesamt.Gesamt`
    :ivar int n_align: number of aligned residues with :py:obj:`~swamp.wrappers.gesamt.Gesamt`

    :example:

    >>> from swamp.wrappers import MapAlign
    >>> mapalign = MapAlign('<workdir>', '<map_a>', '<map_b>')
    >>> mapalign.run()

    """

    def __init__(self, workdir, map_a, map_b, format_a='mapalign', format_b='mapalign', gap_o=-1, gap_e=-0.01,
                 sep_cut=0, pdb_a=None, pdb_b=None, logger=None):

        super(MapAlign, self).__init__(workdir=workdir, logger=logger)

        self.alignment = None
        self.con_sco = None
        self.gap_sco = None
        self.total_sco = None
        self.alignment_length = None
        self.sep_cut = sep_cut
        self.gap_e = gap_e
        self.gap_o = gap_o
        self.qscore = np.nan
        self.rmsd = np.nan
        self.seq_id = np.nan
        self.n_align = np.nan
        self.pdb_a = pdb_a
        self.format_a = format_a
        self.format_b = format_b
        self.map_a = map_a
        self.map_b = map_b
        self.pdb_b = pdb_b

    # ------------------ Some properties ------------------

    @property
    def wrapper_name(self):
        """The name of this `~swamp.wrapper.wrapper.Wrapper` child class (aleigen)"""
        return "map_align"

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
        """Keywords passed through the command line to map_align"""
        return ['-a', self.input_a, '-b', self.input_b, "-sep_cut", str(self.sep_cut), "-gap_e", str(self.gap_e),
                "-gap_o", str(self.gap_o), '-silent']

    @property
    def summary_results(self):
        """A list with the figures of merit obtained"""
        return [self.map_a, self.map_b, self.con_sco, self.gap_sco, self.total_sco, self.alignment_length, self.qscore,
                self.rmsd, self.seq_id, self.n_align]

    @property
    def map_a_id(self):
        """Basename of :py:attr:`~swamp.wrappers.mapaling.MapAlign.map_a` file without extension"""
        return os.path.splitext(os.path.basename(self.map_a))[0]

    @property
    def map_b_id(self):
        """Basename of :py:attr:`~swamp.wrappers.mapaling.MapAlign.map_b` file without extension"""
        return os.path.splitext(os.path.basename(self.map_b))[0]

    @property
    def _reference_mapformat(self):
        """The reference contact map format to be used with the aligner"""
        return 'mapalign'

    # ------------------ Some general methods ------------------

    def _check_alignment_input(self):
        """Check that :py:attr:`~swamp.wrappers.mapaling.MapAlign.map_a` and \
        :py:attr:`~swamp.wrappers.mapaling.MapAlign.map_b` have the \
        :py:attr:`~swamp.wrappers.mapaling.MapAlign._reference_mapformat` and if not convert them using `conkit.io`"""

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
        """Method to extract the figures of merit out the stdout

        :param str logfile: Not in use
        """

        parser = MapAlignParser(logger=self.logger, stdout=self.logcontents)
        parser.parse()

        if parser.error:
            self.error = True
            self.logger.warning('Previous errors while parsing map_align output detected!')
        else:
            self.alignment, self.con_sco, self.gap_sco, self.total_sco, self.alignment_length = parser.summary

    def benchmark(self):
        """Use :py:obj:`~swamp.wrappers.gesamt.Gesamt` to perform a structural alignment between the input structures \
        for benchmarking purposes"""

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
        """Run the :py:attr:`~swamp.wrappers.mapalign.MapAlign.cmd` and store the stdout

        :raises EnvironmentError: if :py:attr:`~swamp.wrappers.wrapper.Wrapper.source` is not found in the system
        """

        if not os.path.isfile(self.source) is None:
            raise EnvironmentError("Cannot find %s executable!" % self.source)

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
