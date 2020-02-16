import os
import swamp
import conkit.io
import numpy as np
from pyjob import cexec
from swamp.parsers import AleigenParser
from swamp.wrappers.mapalign import MapAlign


class AlEigen(MapAlign):
    """Wrapper around al-eigen.

    This class can be used to perform contact map alignments using al-eigen and parse the obtained results. The class
    extends methods and datastructures from `swamp.wrappers.mapalign.MapAlign` and can be used in the same manner.


    :param workdir: working directory
    :type workdir: str
    :param map_a: file name of the map A to be used in the alignment
    :type map_a: str
    :param pdb_a: file name of the pdb A to be used for benchmarking of the alignment
    :type pdb_a: str
    :param format_a: contact prediction file format of map A (default 'aleigen')
    :type format_a: str
    :param map_b: file name of the map B to be used in the alignment
    :type map_b: str
    :param pdb_b: file name of the pdb B to be used for benchmarking of the alignment
    :type pdb_b: str
    :param format_b: contact prediction file format of map B (default 'aleigen')
    :type format_b: str
    :param n_eigen: number of eigen vectors to be used (default '10')
    :type n_eigen: str
    :param get_matrix: if True an alignment matrix is retrieved (default False)
    :type get_matrix: bool
    :param logger: logging interface for the wrapper (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`
    :param eigenvectors: tuple with the file names of both pre-computed eigenvectors used to save time (default None)
    :type eigenvectors: tuple, list
    :ivar c1: number of contacts in map_a
    :ivar _c2: number of contacts in map_b
    :ivar _cmo: contact maximum overlap obtained in this alignment

    :examples

    >>> from swamp.wrappers import AlEigen
    >>> aleigen = AlEigen('<workdir>', '<map_a>', '<map_b>')
    >>> aleigen.run()

    """

    def __init__(self, workdir, map_a, map_b, format_a='aleigen', format_b='aleigen', get_matrix=False, pdb_a=None,
                 pdb_b=None, n_eigen="10", logger=None, eigenvectors=None):

        #if swamp.SRC_ALEIGEN is None:
        #    raise EnvironmentError("Couldn't find aleigen binary!")

        super(AlEigen, self).__init__(workdir=workdir, pdb_a=pdb_a, pdb_b=pdb_b, logger=logger, map_a=map_a,
                                      map_b=map_b, format_a=format_a, format_b=format_b)

        self._c1 = np.nan
        self._c2 = np.nan
        self._cmo = np.nan
        self._get_matrix = get_matrix
        self._n_eigen = n_eigen
        self._eigenvectors = eigenvectors

    # ------------------ Some properties ------------------

    @property
    def eigenvectors(self):
        """Property to store eigenvectors is any is provided"""
        return self._eigenvectors

    @eigenvectors.setter
    def eigenvectors(self, value):
        self._eigenvectors = value

    @property
    def c1(self):
        """Property to store c1"""
        return self._c1

    @c1.setter
    def c1(self, value):
        self._c1 = value

    @property
    def n_eigen(self):
        """Property to store n_eigen"""
        return self._n_eigen

    @n_eigen.setter
    def n_eigen(self, value):
        self._n_eigen = value

    @property
    def c2(self):
        """Property to store c2"""
        return self._c2

    @c2.setter
    def c2(self, value):
        self._c2 = value

    @property
    def cmo(self):
        """Property to store cmo"""
        return self._cmo

    @cmo.setter
    def cmo(self, value):
        self._cmo = value

    @property
    def get_matrix(self):
        """Property to store get_matrix"""
        return self._get_matrix

    @get_matrix.setter
    def get_matrix(self, value):
        self._get_matrix = value

    @property
    def summary_results(self):
        return [self.map_a, self.map_b, self.con_sco, self.c1, self.c2, self.cmo, self.alignment_length, self.qscore,
                self.rmsd, self.seq_id, self.n_align]

    @property
    def _reference_mapformat(self):
        """The reference contact map format to be used with the aligner"""
        return 'aleigen'

    @property
    def alignment_length(self):
        """Length of the contact map alignemnt obtained"""
        if self.alignment is None:
            return None
        else:
            return len(self.alignment.keys())

    @property
    def keywords(self):
        """Keywords used to calculate the alignment"""

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
        return 'aleigen'

    @property
    def source(self):
        """Location of the executable file of al-eigen"""
        return swamp.SRC_ALEIGEN

    # ------------------ Some general methods ------------------

    def get_scores(self, logfile=None):
        """Method to extract the scores out of the stdout

        :param logfile: Not in use
        :returns nothing
        :rtype None
        """

        parser = AleigenParser(logger=self.logger, stdout=self.logcontents)
        parser.parse()

        if parser.error:
            self.error = True
            self.logger.warning('Previous errors while parsing aleigen output detected!')
        else:
            self.alignment, self.con_sco, self.cmo, self.c1, self.c2 = parser.summary

    @staticmethod
    def create_eigen_vectors(cmap, vector_output, map_format='aleigen'):
        """Method to create an eigen vector file from a given contact map

        :param cmap: file name of the input contact map
        :type cmap: str
        :param vector_output: file name to store the output eigen vectors
        :type vector_output: str
        :param map_format: file format of the input contact map file
        :type map_format: str
        :returns nothing
        :rtype None
        """

        if swamp.SRC_WEIGENVECT is None:
            raise EnvironmentError("Cannot find weigenvect binary executable!")

        if map_format != 'aleigen':
            tmpfile = AlEigen.get_tempfile()
            conkit.io.convert(fname_in=cmap, format_in=map_format, fname_out=tmpfile, format_out='aleigen')
            cmd = [swamp.SRC_WEIGENVECT, tmpfile]
        else:
            tmpfile = None
            cmd = [swamp.SRC_WEIGENVECT, cmap]

        with open(vector_output, 'w') as fhandle:
            cexec(cmd, stdout=fhandle)

        if tmpfile is not None:
            os.remove(tmpfile)
