import os
from pyjob import cexec
from swamp.parsers import RefmacParser
from swamp.wrappers.wrapper import Wrapper


class wRefmac(Wrapper):
    """Refmac5 wrapper

    :param str workdir: working directory where refmac will be executed
    :param str pdbin: pdb filename with the placed search model
    :param str mtzin: target's mtz filename
    :param str make_hydr: mset MAKE_HYDR (default: 'N')
    :param str weight_matrix: set WEIGHT MATRIX stdin value for jelly body refinement (default: '0.01')
    :param str ridg_dist_sigm: set RIDG DIST SIGM stdin value for jelly body refinement (default: '0.02')
    :param str ncyc: number of refinement cycles (default: '100')
    :param str phased_mtz: target's mtz filename containing phases (default: None)
    :param bool silent_start: if True, the logger will not display the start banner (default False)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :ivar str rfree: calculated Rfree after refinement
    :ivar str rfactor: calculated Rfactor after refinement
    :ivar tuple rfactor_delta: a tuple with the intial and final Rfactor
    :ivar tuple rfree_delta: a tuple with the intial and final Rfree 
    :ivar tuple bondlenght_delta: a tuple with the intial and final Bond length
    :ivar tuple bondangle_delta: a tuple with the intial and final Bond angle
    :ivar tuple chirvol_delta: a tuple with the intial and final Chirvol
    :ivar str ovarall_CC: overall correlation coefficient between phases at \
    :py:attr:`~swamp.wrappers.wrefmac.wRefmac.phased_mtz` and the refined search model
    :ivar str local_CC: local correlation coefficient between phases at \
    :py:attr:`~swamp.wrappers.wrefmac.wRefmac.phased_mtz` and the refined search model

    :example:

    >>> from swamp.wrappers import wRefmac
    >>> my_refmac = wRefmac('<workdir>', '<pdbin>', '<mtzin>')
    >>> my_refmac.run()
    >>> my_refmac.make_logfile()

    """

    def __init__(self, workdir, pdbin, mtzin, make_hydr='N', weight_matrix="0.01", ncyc="100", ridg_dist_sigm="0.02",
                 logger=None, phased_mtz=None, silent_start=False):

        super(wRefmac, self).__init__(logger=logger, workdir=workdir, silent_start=silent_start)

        self.rfactor = "NA"
        self.rfree = "NA"
        self.rfactor_delta = ("NA", "NA")
        self.rfree_delta = ("NA", "NA")
        self.bondlenght_delta = ("NA", "NA")
        self.bondangle_delta = ("NA", "NA")
        self.chirvol_delta = ("NA", "NA")
        self.local_CC = "NA"
        self.overall_CC = "NA"
        self.make_hydr = make_hydr
        self.ridg_dist_sigm = ridg_dist_sigm
        self.weight_matrix = weight_matrix
        self.ncyc = ncyc
        self.pdbin = pdbin
        self.mtzin = mtzin
        self.phased_mtz = phased_mtz

    @property
    def wrapper_name(self):
        """The name of this wrapper (refmac)"""
        return "refmac"

    @property
    def cmd(self):
        """Command to be executed on the shell"""

        return ["refmac5", 'hklin', self.mtzin, 'hklout', self.mtzout, 'xyzin',
                self.pdbin, 'xyzout', self.pdbout]

    @property
    def summary_results(self):
        """String representation of a summary of the obtained figures of merit"""
        return "Refmac results: Rfactor - %s   Rfree - %s   Local CC - %s   Overall CC - %s" \
               "" % (self.rfactor, self.rfree, self.local_CC, self.overall_CC)

    @property
    def keywords(self):
        """Keywords tos be passed to refmac5 through the stdin"""

        return "RIDG DIST SIGM  %s" % self.ridg_dist_sigm + os.linesep + "MAKE HYDR %s" % \
               self.make_hydr + os.linesep + "WEIGHT MATRIX %s" % self.weight_matrix \
               + os.linesep + "NCYC %s" % self.ncyc + os.linesep + "END"

    def get_scores(self, logfile=None):
        """Parse the logfile and extract scores after refinement using a \
        :py:obj:`~swamp.parsers.refmacparser.RefmacParser` instance
        
        :param None logfile: Not in use (default None)
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
        """Run :py:attr:`~swamp.wrappers.wrefmac.wRefmac.cmd` in the shell"""

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
