import os
from pyjob import cexec
from swamp.parsers import RefmacParser
from swamp.wrappers.wrapper import Wrapper


class wRefmac(Wrapper):
    """Refmac5 wrapper

    :param str workdir: working directory where refmac will be executed
    :param str pdbin: pdb filename with the placed search model
    :param str mtzin: target's mtz filename
    :param str make_hydr: make water molecules during refinement (default: 'N')
    :param str weight_matrix: set WEIGHT MATRIX stdin value for jelly body refinement (default: '0.01')
    :param str ridg_dist_sigm: set RIDG DIST SIGM stdin value for jelly body refinement (default: '0.02')
    :param str ncyc: number of refinement cycles (default: '100')
    :param phased_mtz: target's mtz filename containing phases (default: None)
    :type phased_mtz: str
    :param silent_start: if True, the logger will not display the start banner (default False)
    :type silent_start: bool
    :param logger: logging interface to be used (default None)
    :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
    :ivar rfree: calculated Rfree after refinement
    :type rfree: str
    :ivar rfactor: calculated Rfactor after refinement
    :type rfactor: str

    :example

    >>> from swamp.wrappers import wRefmac
    >>> my_refmac = wRefmac('<workdir>', '<pdbin>', '<mtzin>')
    >>> my_refmac.run()
    >>> my_refmac.make_logfile()

    """

    def __init__(self, workdir, pdbin, mtzin, make_hydr='N', weight_matrix="0.01", ncyc="100", ridg_dist_sigm="0.02",
                 logger=None, phased_mtz=None, silent_start=False):

        super(wRefmac, self).__init__(logger=logger, workdir=workdir, silent_start=silent_start)

        self._rfactor = "NA"
        self._rfree = "NA"
        self._rfactor_delta = ("NA", "NA")
        self._rfree_delta = ("NA", "NA")
        self._bondlenght_delta = ("NA", "NA")
        self._bondangle_delta = ("NA", "NA")
        self._chirvol_delta = ("NA", "NA")
        self._local_CC = "NA"
        self._overall_CC = "NA"
        self._pdbout = os.path.join(self.workdir, "refmac", "refined.pdb")
        self._mtzout = os.path.join(self.workdir, "refmac", "refined.mtz")
        self._make_hydr = make_hydr
        self._ridg_dist_sigm = ridg_dist_sigm
        self._weight_matrix = weight_matrix
        self._ncyc = ncyc
        self._pdbin = pdbin
        self._mtzin = mtzin
        self._phased_mtz = phased_mtz

    @property
    def wrapper_name(self):
        return "refmac"

    @property
    def ridg_dist_sigm(self):
        return self._ridg_dist_sigm

    @ridg_dist_sigm.setter
    def ridg_dist_sigm(self, value):
        self._ridg_dist_sigm = value

    @property
    def weight_matrix(self):
        return self._weight_matrix

    @weight_matrix.setter
    def weight_matrix(self, value):
        self._weight_matrix = value

    @property
    def ncyc(self):
        return self._ncyc

    @ncyc.setter
    def ncyc(self, value):
        self._ncyc = value

    @property
    def mtzin(self):
        return self._mtzin

    @mtzin.setter
    def mtzin(self, value):
        self._mtzin = value

    @property
    def phased_mtz(self):
        return self._phased_mtz

    @phased_mtz.setter
    def phased_mtz(self, value):
        self._phased_mtz = value

    @property
    def pdbin(self):
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def make_hydr(self):
        return self._make_hydr

    @make_hydr.setter
    def make_hydr(self, value):
        self._make_hydr = value

    @property
    def rfactor(self):
        return self._rfactor

    @rfactor.setter
    def rfactor(self, value):
        self._rfactor = value

    @property
    def rfactor_delta(self):
        return self._rfactor_delta

    @rfactor_delta.setter
    def rfactor_delta(self, value):
        self._rfactor_delta = value

    @property
    def rfree(self):
        return self._rfree

    @rfree.setter
    def rfree(self, value):
        self._rfree = value

    @property
    def rfree_delta(self):
        return self._rfree_delta

    @rfree_delta.setter
    def rfree_delta(self, value):
        self._rfree_delta = value

    @property
    def bondlenght_delta(self):
        return self._bondlenght_delta

    @bondlenght_delta.setter
    def bondlenght_delta(self, value):
        self._bondlenght_delta = value

    @property
    def bondangle_delta(self):
        return self._bbondangle_delta

    @bondangle_delta.setter
    def bondangle_delta(self, value):
        self._bondangle_delta = value

    @property
    def chirvol_delta(self):
        return self._chirvol_delta

    @chirvol_delta.setter
    def chirvol_delta(self, value):
        self._chirvol_delta = value

    @property
    def local_CC(self):
        """Local CC with between the mtzfile and the refined pdbout"""
        return self._local_CC

    @local_CC.setter
    def local_CC(self, value):
        self._local_CC = value

    @property
    def overall_CC(self):
        """Overall CC between the mtzfile and the refined pdbout"""
        return self._overall_CC

    @overall_CC.setter
    def overall_CC(self, value):
        self._overall_CC = value

    @property
    def cmd(self):
        """Command to be executed on the shell"""

        return ["refmac5", 'hklin', self.mtzin, 'hklout', self.mtzout, 'xyzin',
                self.pdbin, 'xyzout', self.pdbout]

    @property
    def summary_results(self):
        """String with a summary of the obtained figures of merit"""
        return "Refmac results: Rfactor - %s   Rfree - %s   Local CC - %s   Overall CC - %s" \
               "" % (self.rfactor, self.rfree, self.local_CC, self.overall_CC)

    @property
    def keywords(self):
        """Keywords tos be passed to refmac5 through the stdin"""

        return "RIDG DIST SIGM  %s" % self.ridg_dist_sigm + os.linesep + "MAKE HYDR %s" % \
               self.make_hydr + os.linesep + "WEIGHT MATRIX %s" % self.weight_matrix \
               + os.linesep + "NCYC %s" % self.ncyc + os.linesep + "END"

    def get_scores(self, logfile=None):
        """Parse the logfile and extract scores after refinement
        
        :param logfile: Not in use (default None)
        :type logfile: None
        :returns nothing
        :rtype None
        """

        parser = RefmacParser(self.logcontents, logger=self.logger)
        parser.parse()

        # If there was an error
        if parser.error:
            self.error = True
            self.logger.warning('Previous detected while parsing refmac output!')
            return

        self.rfactor, self.rfree, self.rfactor_delta, self.rfree_delta, self.bondlenght_delta, \
        self.bondangle_delta, self.chirvol_delta = parser.summary

        # If everything is ok so far, get the CC
        if self.phased_mtz is not None:
            self.local_CC, self.overall_CC = self.get_cc(os.path.join(self.workdir, "phenix"), self.phased_mtz,
                                                         self.pdbout)

    def run(self):
        """Run refmac5 with specified parameters"""

        if not self.silent_start:
            self.logger.info(self.wrapper_header)
        self.logger.info('Running %s cycles of refinement' % self.ncyc)

        self.make_workdir()
        tmp_dir = os.getcwd()
        os.chdir(self.workdir)

        self.logger.debug(" ".join(self.cmd))
        self.logcontents = cexec(self.cmd, stdin=self.keywords, permit_nonzero=True)

        self.get_scores()

        if not self.error and (not os.path.isfile(self.pdbout) or not os.path.isfile(self.mtzout)):
            self.logger.error("refmac output pdb/mtz missing!")
            self.error = True

        self.logger.info(self.summary_results)

        os.chdir(tmp_dir)
