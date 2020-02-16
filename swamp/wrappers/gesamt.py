import os
import gemmi
import itertools
import threading
import numpy as np
import pandas as pd
from pyjob import cexec
from swamp.wrappers import Wrapper
from swamp.utils import ThreadResults, invert_hiearchy, get_tempfile


class Gesamt(Wrapper):
    """Wrapper around gesamt

    :param workdir: working directory
    :type workdir: str, None
    :param mode: specify the type of task to be done with gesamt
    :type mode: str
    :param pdbin: input pdb file name. If an iterable is given, assumes multiple structure alignment
    :type pdbin: str, tuple, list
    :param gesamt_archive: location of the gesamt archive to be scanned (default None)
    :type gesamt_archive: str, None
    :param pdbout: output pdb file name (default None)
    :type pdbout: str, None
    :param nthreads: number of threads to be used by gesamt (default 1)
    :type nthreads: int
    :param hits_out: file name of the .hits output file (default None)
    :type hits_out: str, None
    :param min1: argument to be passed to gesamt as -min1 (default 0.1)
    :type min1: float
    :param min2: argument to be passed to gesamt as -min2 (default 0.1)
    :type min2: float
    :param pdb_archive: location of the pdb archive to be used in gesamt (default None)
    :type pdb_archive: str, None
    :param logger: logging interface for the wrapper
    :type logger: :object:`swamp.logger.swamplogger.SwampLogger`
    :ivar error: if True an error has occurred along the process
    :type error: bool
    """

    def __init__(self, workdir, mode, pdbin, gesamt_archive=None, pdbout=None, nthreads=1, hits_out=None, min1=0.1,
                 min2=0.1, pdb_archive=None, logger=None):

        super(Gesamt, self).__init__(workdir=workdir, logger=logger)

        self._summary_results = None
        self._qscore = None
        self._rmsd = None
        self._seq_id = None
        self._n_align = None
        self._mode = mode
        self._pdbin = pdbin
        self._gesamt_archive = gesamt_archive
        self._pdbout = pdbout
        self._nthreads = nthreads
        self._hits_out = hits_out
        self._pdb_archive = pdb_archive
        self._min1 = min1
        self._min2 = min2

    # ------------------ Properties ------------------

    @property
    def wrapper_name(self):
        return "gesamt"

    @property
    def pdbin(self):
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def min2(self):
        return self._min2

    @min2.setter
    def min2(self, value):
        self._min2 = value

    @property
    def hits_out(self):
        return self._hits_out

    @hits_out.setter
    def hits_out(self, value):
        self._hits_out = value

    @property
    def min1(self):
        return self._min1

    @min1.setter
    def min1(self, value):
        self._min1 = value

    @property
    def nthreads(self):
        return self._nthreads

    @nthreads.setter
    def nthreads(self, value):
        self._nthreads = value

    @property
    def pdbout(self):
        return self._pdbout

    @pdbout.setter
    def pdbout(self, value):
        self._pdbout = value

    @property
    def gesamt_archive(self):
        return self._gesamt_archive

    @gesamt_archive.setter
    def gesamt_archive(self, value):
        self._gesamt_archive = value

    @property
    def pdb_archive(self):
        return self._pdb_archive

    @pdb_archive.setter
    def pdb_archive(self, value):
        self._pdb_archive = value

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        self._mode = value

    @property
    def qscore(self):
        """Qscore obtained through the structural alignment"""
        return self._qscore

    @qscore.setter
    def qscore(self, value):
        self._qscore = value

    @property
    def n_align(self):
        """Number of aligned residues as obtained through the structural alignment"""
        return self._n_align

    @n_align.setter
    def n_align(self, value):
        self._n_align = value

    @property
    def seq_id(self):
        """Sequence identity as obtained through the structural alignment"""
        return self._seq_id

    @seq_id.setter
    def seq_id(self, value):
        self._seq_id = value

    @property
    def rmsd(self):
        """RMSD obtained through the structural alignment"""
        return self._rmsd

    @rmsd.setter
    def rmsd(self, value):
        self._rmsd = value

    @property
    def summary_results(self):
        """A summary of the results obtained through the structural alignment"""
        return self._summary_results

    @summary_results.setter
    def summary_results(self, value):
        self._summary_results = value

    @property
    def keywords(self):
        """No keywords are used in gesamt through stdin"""
        return None

    @property
    def _SRC_GESAMT(self):
        """Location of the executable binary gesatm file"""
        return os.path.join(os.environ['CCP4'], 'bin', 'gesamt')

    @property
    def cmd(self):
        """Command to be executed in the shell"""

        # Create an archive
        if self.mode == "make-archive":
            if self.gesamt_archive is None or self.pdb_archive is None:
                self.logger.error('Impossible to create archive, need to provide location!')
                self.error = True
                return None
            return [self._SRC_GESAMT, '--make-archive', self.gesamt_archive, '-pdb', self.pdb_archive]

        # Align structures
        elif self.mode == "alignment":
            cmd = [self._SRC_GESAMT]
            if not isinstance(self.pdbin, list) and not isinstance(self.pdbin, tuple):
                self.pdbin = list(self.pdbin)
            if len(self.pdbin) < 2:
                self.error = True
                self.logger.error('Need to provide at least two structures for alignment!')
                return None
            for fname in self.pdbin:
                cmd.append(fname)
            if self.pdbout is not None:
                cmd.append('-o')
                cmd.append(self.pdbout)
            return cmd

        # Scan an archive
        elif self.mode == "search-archive":
            cmd = [self._SRC_GESAMT, self.pdbin, '-archive', self.gesamt_archive, '-o', self.hits_out]
            self._filthy_files.append(self.hits_out)
            if self.nthreads is not None:
                cmd.append('-nthreads=%s' % self.nthreads)
            if self.min1 is not None:
                cmd.append('-min1=%s' % self.min1)
            if self.min2 is not None:
                cmd.append('-min2=%s' % self.min2)
            return cmd

        # Unkown mode
        else:
            self.error = True
            self.logger.error("Unkown mode! %s" % self.mode)
            self.logger.warning('Allowed modes are: %s' % ",".join(['alignment', "search-archive", "make-archive"]))
            return None

    # ------------------ General methods ------------------

    def get_scores(self, logfile=None):
        """Method to extract the scores and figures of merit out of the stdout

        :param logfile: No in use
        :type logfile: None
        :returns nothing
        :rtype None
        """

        if self.mode == "search-archive":
            self.summary_results = self.parse_hitfile(self.hits_out)
        elif self.mode == "alignment":
            self.qscore, self.rmsd, self.seq_id, self.n_align = self.parse_stdout(self.logcontents, len(self.pdbin))

    def run(self):
        """Run gesamt and parse the results

        :raises TypeError: command line command is malformed, shouldn't occur please report if this happens
        :returns nothing
        :rtype None
        """

        if self.mode == 'make-archive' and not os.path.isfile(self.gesamt_archive):
            os.makedirs(self.gesamt_archive)

        try:
            self.logger.debug(" ".join(self.cmd))
        except TypeError:
            self.logger.error('CMD IS NOT ITERABLE, THIS SHOULD NOT HAPPEN PLEASE REPORT ISSUE')
            self.logger.debug(self.cmd)
            self.logger.debug(self.pdbin)
            self.logger.debug(self.mode)
            self.error = True
            return
        self.logcontents = cexec(self.cmd, permit_nonzero=True)

        if self.logcontents == b'':
            self.error = True
            self.logger.error("Something went wrong, no gesamt stdout! Exiting now...")
            return

        elif 'DISSIMILAR' in self.logcontents:
            self.error = True
            self.logger.warning("%s are to dissimilar to be aligned!" % (" ".join(self.pdbin)))
            return

        elif 'ALIGNMENT ERROR 2' in self.logcontents:
            self.error = True
            self.logger.warning("%s alignment returned ALIGNMENT ERROR 2!" % (" ".join(self.pdbin)))
            return

        self.get_scores()
        self._cleanup_files()

    # ------------------ Static methods ------------------

    @staticmethod
    def parse_stdout(stdout, n_models):
        """Method to retrieve qscore, rmsd, sequence identity and no. of aligned residues from gesamt stdout

        :param stdout: gesamt stdout for a given command
        :type stdout: str
        :param n_models: number of models that were used in the structural alignment
        :type n_models: int
        :returns qscore, rmsd, sequence identity and no. of aligned residues
        :rtype tuple
        """

        if n_models == 2:
            qscore_mark = "Q-score"
            rmsd_mark = "RMSD"
            n_align_mark = "Aligned residues"
            seqid_mark = "Sequence Id"
        else:
            qscore_mark = "quality Q"
            rmsd_mark = "r.m.s.d"
            n_align_mark = "Nalign"
            seqid_mark = "this will not be here"

        qscore = np.nan
        rmsd = np.nan
        n_align = np.nan
        seq_id = np.nan
        for line in stdout.split("\n"):
            if len(line.split()) != 0 and line.split()[0] != "#":
                if qscore_mark in line and qscore is np.nan:
                    qscore = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif rmsd_mark in line and rmsd is np.nan:
                    rmsd = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif n_align_mark in line and n_align is np.nan:
                    n_align = int(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
                elif seqid_mark in line and seq_id is np.nan:
                    seq_id = float(line.rstrip().lstrip().split(":")[-1].split()[0].rstrip().lstrip())
        return qscore, rmsd, seq_id, n_align

    @staticmethod
    def parse_hitfile(fname):
        """Method to parse a gesamt .hit output file

        :param fname: file name of the .hit output file
        :type fname: str
        :returns df: a dataframe with the results contained in the hit file
        :rtype :object:`pandas.DataFrame`
        """

        df = []
        with open(fname, "r") as fhandle:
            for line in fhandle:
                if line[0] != "#":
                    line = line[20:].split()
                    df.append([line[-6], line[-5], line[-4], line[-3], line[-2], line[-1]])
        df = pd.DataFrame(df)
        df.columns = ["qscore", "rmsd", "seq_id", "n_align", "n_res", "fname"]
        return df

    @staticmethod
    def get_pairwise_qscores(stdout):
        """Method to get the pairwise qscores of a given alignmnet between several structures in an ensemble

        :param stdout: gesamt stdout for the command
        :type stdout: str
        :returns: qscores_dict: a dictionary with the pairwise qscores for each of the structures in the alignment
        :rtype dict
        """

        qscores_dict = {}
        structure_id_dict = {}
        qscores_mark = "(o) pairwise Q-scores"
        file_mark = "... reading file"
        rmsd_mark = "(o) pairwise r.m.s.d."
        is_qscores = False

        for line in stdout.split("\n"):

            # Store file names and structure ids
            if file_mark in line:
                fname = line.split("'")[1]
                structure_id = "S%s" % str(len(qscores_dict.keys()) + 1).zfill(3)
                qscores_dict[fname] = None
                structure_id_dict[structure_id] = fname
            # Qscores will start appearing now
            elif qscores_mark in line:
                is_qscores = True
            # Store the qscore in the dictionary
            elif is_qscores and line.split("|")[0].rstrip().lstrip() in structure_id_dict.keys():
                structure_id = line.split("|")[0].rstrip().lstrip()
                idx = int(structure_id[1:])
                qscores_dict[structure_id_dict[structure_id]] = float(line.split()[idx].rstrip().lstrip())
            # If we reach the rmsd mark, break the loop
            elif rmsd_mark in line:
                break

        return qscores_dict

    @staticmethod
    def get_optimum_alignment(pdbfiles, nthreads=1, logger=None):
        """Method to get the optimum alignment between a set given pdb files.

        This method considers all the structure arrangements: it screens for all the possible sets of reference
        structures and inverted structures possible). This can be computationally demanding so multi-threading
        is possible.

        :param pdbfiles: the set of pdb files to be aligned
        :type pdbfiles: tuple, list
        :param nthreads: number of threads to be used when screening all possible arrangements (default 1)
        :type nthreads: int
        :param logger: logging interface to be used (default None)
        :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
        :returns the gesamt instance for the optimal alignment and a pdb hierarchy with the aligned structures ensembled
        :rtype tuple
        """

        def _align(idx, combination, logger):
            """Align the specified combination of pdbfiles. It uses a semaphore for multiple thread coordination.

            :param idx: index of the combination to be aligned (used as identifier)
            :type idx: int
            :param combination: the combination of files to be aligned
            :type combination: tuple
            :returns nothing
            :rtype None
            """

            semaphore.acquire()
            gesamt = Gesamt(workdir=None, logger=logger, pdbin=[fname for fname in combination], mode="alignment")
            gesamt.run()
            if gesamt.qscore:
                results.register((idx, gesamt.qscore))
            else:
                results.register((idx, 0.0))
            semaphore.release()

        if len(pdbfiles) == 1:
            return None, gemmi.read_structure(pdbfiles[0])

        # Initiate values
        input_list = []
        tmp_pdbout = get_tempfile()
        tmp_files = [tmp_pdbout]
        child_threads = []
        semaphore = threading.Semaphore(value=nthreads)
        results = ThreadResults()

        # Write the inverted pdb files
        for pdbfile in pdbfiles:
            hierarchy = gemmi.read_structure(pdbfile)
            inverted = invert_hiearchy(hierarchy)
            tmp_inverted_file = get_tempfile()
            inverted.write_minimal_pdb(tmp_inverted_file)
            tmp_files.append(tmp_inverted_file)
            input_list.append((pdbfile, tmp_inverted_file))

        # Determine the alignment with the highest qscore
        all_combinations = list(itertools.product(*input_list))
        for idx, combination in enumerate(all_combinations):
            t = threading.Thread(target=_align, args=(idx, combination, logger,))
            child_threads.append(t.getName())
            t.start()
        mainthread = threading.current_thread()
        for t in threading.enumerate():
            if t is not mainthread and t.getName() in child_threads:
                t.join()

        ranked_alignemnts = sorted(results.value, key=lambda x: x[1], reverse=True)
        optimal_alignment = Gesamt(workdir=None, logger=logger, pdbin=all_combinations[ranked_alignemnts[0][0]],
                                   mode="alignment", pdbout=tmp_pdbout)
        optimal_alignment.run()
        if not optimal_alignment.error:
            if not os.path.isfile(tmp_pdbout):
                logger.debug('Not found! %s' % tmp_pdbout)
                logger.debug(all_combinations[ranked_alignemnts[0][0]])
                logger.debug(os.listdir(os.environ['CCP4_SCR']))
            optimal_ensemble = gemmi.read_structure(tmp_pdbout)
        else:
            optimal_ensemble = None

        # Remove tmp files and return optimal alignment
        for fname in tmp_files:
            if os.path.isfile(fname):
                os.remove(fname)

        return optimal_alignment, optimal_ensemble
