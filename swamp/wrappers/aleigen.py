import os
import swamp
import conkit.io
import numpy as np
from pyjob import cexec
from swamp.utils import get_tempfile
from swamp.parsers import AleigenParser
from swamp.wrappers.mapalign import MapAlign


class AlEigen(MapAlign):
    """Wrapper around al-eigen.

    This class can be used to perform contact map alignments using al-eigen and parse the obtained results. The class \
    extends methods and datastructures from :py:obj:`~swamp.wrappers.mapalign.MapAlign` and can be used in the same \
    manner.

    :param str workdir: working directory for this :py:obj:`~swamp.wrappers.aleigen.AlEigen` instance
    :param str map_a: file name of the map A to be used in the alignment
    :param str pdb_a: file name of the pdb A to be used for benchmarking of the alignment
    :param str format_a: contact prediction file format of map A (default 'aleigen')
    :param str map_b: file name of the map B to be used in the alignment
    :param str pdb_b: file name of the pdb B to be used for benchmarking of the alignment
    :param str format_b: contact prediction file format of map B (default 'aleigen')
    :param str n_eigen: number of eigen vectors to be used in the alignment (default '10')
    :param bool get_matrix: if True an alignment matrix is retrieved (default False)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :param tuple eigenvectors: tuple with the file names of both pre-computed eigenvectors (default None)
    :ivar str c1: number of contacts in map_a
    :ivar str c2: number of contacts in map_b
    :ivar str cmo: contact maximum overlap obtained in this alignment

    :example:

    >>> from swamp.wrappers import AlEigen
    >>> aleigen = AlEigen('<workdir>', '<map_a>', '<map_b>')
    >>> aleigen.run()

    """

    def __init__(self, workdir, map_a, map_b, format_a='aleigen', format_b='aleigen', get_matrix=False, pdb_a=None,
                 pdb_b=None, n_eigen="10", logger=None, eigenvectors=None):

        if swamp.SRC_ALEIGEN is None:
            raise EnvironmentError("Couldn't find aleigen binary!")

        super(AlEigen, self).__init__(workdir=workdir, pdb_a=pdb_a, pdb_b=pdb_b, logger=logger, map_a=map_a,
                                      map_b=map_b, format_a=format_a, format_b=format_b)

        self.c1 = np.nan
        self.c2 = np.nan
        self.cmo = np.nan
        self.get_matrix = get_matrix
        self.n_eigen = n_eigen
        self.eigenvectors = eigenvectors

    # ------------------ Some properties ------------------

    @property
    def summary_results(self):
        """A list with a summary of the obtained figures of merit"""
        return [self.map_a, self.map_b, self.con_sco, self.c1, self.c2, self.cmo, self.alignment_length, self.qscore,
                self.rmsd, self.seq_id, self.n_align]

    @property
    def _reference_mapformat(self):
        """The reference contact map format to be used with the aligner"""
        return 'aleigen'

    @property
    def keywords(self):
        """Keywords passed through the command line to aleigen"""

        keywords = [self.input_a, self.input_b, self.n_eigen]

        if self.eigenvectors is not None:
            keywords.append('-e')
            keywords.append(self.eigenvectors[0])
            keywords.append(self.eigenvectors[1])
        if self.get_matrix:
            keywords.append('-m')

        return keywords

    @property
    def wrapper_name(self):
        """The name of this wrapper (aleigen)"""
        return 'aleigen'

    @property
    def source(self):
        """Location of the executable file of aleigen"""
        return swamp.SRC_ALEIGEN

    # ------------------ Some general methods ------------------

    def get_scores(self, logfile=None):
        """Method to extract the figures of merit out the stdout

        :param logfile: Not in use
        """

        parser = AleigenParser(logger=self.logger, stdout=self.logcontents)
        parser.parse()

        if parser.error:
            self.error = True
            self.logger.warning('Previous errors while parsing aleigen output detected!')
        else:
            self.alignment, self.alignment_length, self.con_sco, self.cmo, self.c1, self.c2 = parser.summary

    @staticmethod
    def create_eigen_vectors(cmap, vector_output, map_format='aleigen'):
        """Create an eigen vector file from a given contact map

        :param str cmap: file name of the input contact map
        :param str vector_output: file name to store the output eigen vectors
        :param str map_format: file format of the input contact map file (default 'aleigen')
        """

        if swamp.SRC_WEIGENVECT is None:
            raise EnvironmentError("Cannot find weigenvect binary executable!")

        if map_format != 'aleigen':
            tmpfile = get_tempfile()
            conkit.io.convert(fname_in=cmap, format_in=map_format, fname_out=tmpfile, format_out='aleigen')
            cmd = [swamp.SRC_WEIGENVECT, tmpfile]
        else:
            tmpfile = None
            cmd = [swamp.SRC_WEIGENVECT, cmap]

        with open(vector_output, 'w') as fhandle:
            cexec(cmd, stdout=fhandle)

        if tmpfile is not None:
            os.remove(tmpfile)

