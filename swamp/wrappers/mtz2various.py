import os
from pyjob import cexec
from swamp.parsers import MtzParser
from swamp.wrappers.wrapper import Wrapper


class Mtz2Various(Wrapper):
    """Wrapper around mtz2various

    :param str workdir: working directory
    :param str mtzin: mtz file name
    :param str hklout: hkl output file name
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :param bool silent_start: if True, the logger will not display the start banner (default False)
    :ivar bool error: if True an error has occurred along the process

    :example:

    >>> from swamp.wrappers import Mtz2Various
    >>> my_mtz2various = Mtz2Various('<workdir>' '<mtzin>', '<hklout>')
    >>> my_mtz2various.run()
    >>> my_mtz2various.make_logfile()
    """

    def __init__(self, workdir, mtzin, hklout, logger=None, silent_start=False):

        super(Mtz2Various, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self.mtzin = mtzin
        self.hklout = hklout

    @property
    def keywords(self):
        """Keywords to pass to mtz2various through stdin"""

        mtz_head = MtzParser(self.mtzin)
        mtz_head.parse()

        # We prefer I over F
        if mtz_head.i is not None and mtz_head.sigi is not None:
            return 'LABIN I=%s SIGI=%s FREE=%s' % (
                mtz_head.i, mtz_head.sigi, mtz_head.free) + os.linesep + "OUTPUT SHELX" + os.linesep + "END"
        else:
            return 'LABIN FP=%s SIGFP=%s FREE=%s' % (mtz_head.f, mtz_head.sigf, mtz_head.free) \
                   + os.linesep + "OUTPUT SHELX" + os.linesep + "FSQUARED" + os.linesep + "END"


    @property
    def summary_results(self):
        """No figures of merit are obtained with mtz2various"""
        return None

    @property
    def wrapper_name(self):
        """The name of this `~swamp.wrapper.wrapper.Wrapper` child class (mtz2various)"""
        return "mtz2various"

    @property
    def cmd(self):
        """Command to be executed on the shell"""
        return [self.source, 'HKLIN', self.mtzin, 'HKLOUT', self.hklout]

    def get_scores(self, logfile=None):
        """Abstract method to get scores (not implemented in this class)"""
        pass

    def _run(self):
        """Run the :py:attr:`~swamp.wrappers.mtz2various.Mtz2Various.cmd` and store the stdout"""

        self.logger.info(self.wrapper_header)
        self.logger.debug(" ".join(self.cmd))
        self.make_workdir()
        self.logcontents = cexec(self.cmd, stdin=self.keywords)
