import os
from pyjob import cexec
from swamp.parsers import MtzParser
from swamp.wrappers import Wrapper


class Mtz2Various(Wrapper):
    """Wrapper around mtz2various

     :param workdir: working directory
     :type workdir: str, None
     :param mtzin: mtz file name
     :type mtzin: str
     :param hklout: hkl output file name
     :type hklout: str
     :param logger: logging interface to be used (default None)
     :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
     :param silent_start: if True, the logger will not display the start banner (default False)
     :type silent_start: bool
     :ivar error: if True an error has occurred along the process
     :type error: bool

    :example

    >>> from swamp.wrappers import Mtz2Various
    >>> my_mtz2various = Mtz2Various('<workdir>' '<mtzin>', '<hklout>')
    >>> my_mtz2various.run()
    >>> my_mtz2various.make_logfile()
    """

    def __init__(self, workdir, mtzin, hklout, logger=None, silent_start=False):

        super(Mtz2Various, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self._mtzin = mtzin
        self._hklout = hklout

    @property
    def mtzin(self):
        return self._mtzin

    @mtzin.setter
    def mtzin(self, value):
        self._mtzin = value

    @property
    def hklout(self):
        return self._hklout

    @hklout.setter
    def hklout(self, value):
        self._hklout = value

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
    def wrapper_name(self):
        return "mtz2various"

    @property
    def cmd(self):
        """Command to be executed on the shell"""
        return ["mtz2various", 'HKLIN', self.mtzin, 'HKLOUT', self.hklout]

    def get_scores(self, logfile=None):
        """Abstract method to get scores (not implemented in this class)"""
        pass

    def run(self):
        """Run the mtz2various and store the stdout"""

        self.logger.info(self.wrapper_header)
        self.logger.debug(" ".join(self.cmd))
        self.make_workdir()
        self.logcontents = cexec(self.cmd, stdin=self.keywords)
