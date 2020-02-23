import os
import gemmi
import itertools
import threading
from pyjob import cexec
from swamp.parsers import GesamtParser
from swamp.wrappers.wrapper import Wrapper
from swamp.utils import ThreadResults, invert_hiearchy, get_tempfile


class Gesamt(Wrapper):
    """Wrapper around gesamt

    This class can be used to perform structural alignments using gesamt and parse the obtained results.

    :param str workdir: working directory :py:obj:`~swamp.wrappers.gesamt.Gesamt` instance
    :type workdir: str, None
    :param str mode: specify the type of task to be done with gesamt
    :param pdbin: input pdb file name. If an iterable is given, assumes multiple structure alignment
    :type pdbin: str, tuple, list
    :param gesamt_archive: location of the gesamt archive to be scanned (default None)
    :type gesamt_archive: str, None
    :param pdbout: output pdb file name (default None)
    :type pdbout: str, None
    :param int nthreads: number of threads to be used by gesamt (default 1)
    :param hits_out: file name of the .hits output file (default None)
    :type hits_out: str, None
    :param float min1: argument to be passed to gesamt as -min1 (default 0.1)
    :param float min2: argument to be passed to gesamt as -min2 (default 0.1)
    :param pdb_archive: location of the pdb archive to be used in gesamt (default None)
    :type pdb_archive: str, None
    :param `swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper
    :ivar bool error: if True an error has occurred along the process
    :ivar float qscore: qscore as reported by gesamt
    :ivar float rmsd: the obtained rmsd as reported by gesamt
    :ivar float seq_id: sequence identity between the input structures
    :ivar int n_align: number of aligned residues
    """

    def __init__(self, workdir, mode, pdbin=None, gesamt_archive=None, pdbout=None, nthreads=1, hits_out=None, min1=0.1,
                 min2=0.1, pdb_archive=None, logger=None):

        super(Gesamt, self).__init__(workdir=workdir, logger=logger)

        self.summary_results = None
        self.qscore = None
        self.rmsd = None
        self.seq_id = None
        self.n_align = None
        self.mode = mode
        self.pdbin = pdbin
        self.gesamt_archive = gesamt_archive
        self._pdbout = pdbout
        self.nthreads = nthreads
        self.hits_out = hits_out
        self.pdb_archive = pdb_archive
        self.min1 = min1
        self.min2 = min2

        if self.pdbin is None and self.mode != 'make-archive':
            self.logger.error('Need to provide at least one PDB file for gesamt!')
            self.error = True
        elif self.mode not in self.allowed_modes:
            self.logger.error(
                'Selected mode was not recognised! Allowed modes are: {}'.format(' '.join(self.allowed_modes)))
            self.error = True

    # ------------------ Properties ------------------

    @property
    def pdbout(self):
        """Setter for :py:attr:`~swamp.wrapper.gesamt.Gesamt.pdbout`"""
        return self._pdbout

    @pdbout.setter
    def pdbout(self, value):
        self._pdbout = value

    @property
    def allowed_modes(self):
        """A list with the allowed values for :py:attr:`~swmap.wrappers.gesamt.Gesamt.mode`"""
        return ['make-archive', 'alignment', 'search-archive']

    @property
    def wrapper_name(self):
        """The name of this `~swamp.wrapper.wrapper.Wrapper` child class (gesamt)"""
        return "gesamt"

    @property
    def summary_results(self):
        """A summary of the figures of merit obtained through the structural alignment"""
        return self._summary_results

    @summary_results.setter
    def summary_results(self, value):
        self._summary_results = value

    @property
    def keywords(self):
        """No keywords are used in gesamt through stdin"""
        return None

    @property
    def cmd(self):
        """Command to be executed on :py:func:`~swamp.wrappers.wrapper.Wrapper.run`"""

        # Create an archive
        if self.mode == "make-archive":
            if self.gesamt_archive is None or self.pdb_archive is None:
                self.logger.error('Impossible to create archive, need to provide location!')
                self.error = True
                return None
            return [self.source, '--make-archive', self.gesamt_archive, '-pdb', self.pdb_archive]

        # Align structures
        elif self.mode == "alignment":
            cmd = [self.source]
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
        else:
            cmd = [self.source, self.pdbin, '-archive', self.gesamt_archive, '-o', self.hits_out]
            self._filthy_files.append(self.hits_out)
            if self.nthreads is not None:
                cmd.append('-nthreads=%s' % self.nthreads)
            if self.min1 is not None:
                cmd.append('-min1=%s' % self.min1)
            if self.min2 is not None:
                cmd.append('-min2=%s' % self.min2)
            return cmd

    # ------------------ General methods ------------------

    def get_scores(self, logfile=None):
        """Use a :py:obj:`~swamp.parsers.gesamtparser.GesamtParser` instance to extract the scores and figures of \
        merit out of the :py:attr:`~swamp.wrappers.wrapper.Wrapper.logcontents` and \
        :py:attr:`~swamp.wrappers.gesamt.Gesamt.hits_out`

        :param logfile: Not in use
        """

        parser = GesamtParser(fname=self.hits_out, mode=self.mode, stdout=self.logcontents)
        parser.parse()

        if self.mode == "search-archive":
            self.summary_results = parser.summary
        elif self.mode == "alignment":
            self.qscore, self.rmsd, self.seq_id, self.n_align = parser.summary

    def _run(self):
        """Run :py:attr:`~swamp.wrappers.gesamt.Gesamt.cmd` and parse the results using \
        :py:func:`~swamp.wrappers.gesamt.Gesamt.get_scores`

        :raises TypeError: command line command is malformed, please report if this occurs
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
    def get_optimum_alignment(pdbfiles, nthreads=1, logger=None):
        """Method to get the optimum alignment between a set given pdb files.

        This method considers all the structure arrangements: it screens for all the possible sets of reference
        structures and inverted structures possible). This can be computationally demanding so multi-threading
        is possible.

        :param pdbfiles: the set of pdb files to be aligned
        :type pdbfiles: tuple, list
        :param int nthreads: number of threads to be used when screening all possible arrangements (default 1)
        :param logger: logging interface to be used (default None)
        :type logger: None, :py:obj:`~swamp.logger.swamplogger.SwampLogger`
        :returns: the :py:obj:`~swamp.wrappers.gesamt.Gesamt` instance instance for the optimal alignment and a \
         :py:obj:`gemmi.Structure` hierarchy with the aligned structuresas an ensemble (tuple)
        """

        def _align(idx, combination, logger):
            """Create a :py:obj:`~swamp.wrappers.gesamt.Gesamt` instance and align the specified combination of pdb \
            files. It uses a :py:obj:`threading.semaphore` for multiple thread coordination.

            :param int idx: index of the combination to be aligned (used as identifier)
            :param tuple combination: the combination of files to be aligned
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
