import os
import math
from pyjob import cexec
from shutil import copyfile
from swamp.wrappers.wrapper import Wrapper
from swamp.parsers.shelxeparser import ShelxeParser
from swamp.wrappers.mtz2various import Mtz2Various


class Shelxe(Wrapper):
    """Wrapper around shelxe

     :param workdir: working directory
     :type workdir: str, None
     :param mtzin: target structure mtz file name
     :type mtzin: str
     :param pdbin: pdb file name with the refined search model
     :type pdbin: str
     :param solvent: estimated solvent content for the crystal (default 0.45)
     :type solvent: str
     :param autotracing_ncyc: number of autotracing cycles to run
     :type autotracing_ncyc: str
     :param density_modif_ncyc: number of density modification cycles to run
     :type density_modif_ncyc: str
     :param use_f: if True use F labels in the mtz file
     :type use_f: bool
     :param prune_residues: if True set -o option (default True)
     :type prune_residues: bool
     :param apply_ncs: sets parameter -n (default True)
     :type apply_ncs: bool
     :param init_time: parameter to set -t (default True)
     :type init_time: str
     :param nreflections: number of reflections in the mtz file (default 40000)
     :type nreflections: int
     :param alphahelix: if True sets the -t parameter (default True)
     :type alphahelix: bool
     :param resolution: resolution of the data in the mtz file (default 2.0)
     :type resolution: float
     :param silent_start: if True, the logger will not display the start banner (default False)
     :type silent_start: bool
     :param logger: logging interface to be used (default None)
     :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
     :ivar error: if True an error has occurred along the process
     :type error: bool

    :examples

    >>> from swamp.wrappers.shelxe import Shelxe
    >>> my_shelxe = Shelxe('<workdir>', '<mtzin>', '<pdbin>')
    >>> my_shelxe.run()
    >>> my_shelxe.make_logfile()

    """

    def __init__(self, workdir, pdbin, mtzin, solvent='0.45', autotracing_ncyc='20', density_modif_ncyc='10',
                 use_f=True, prune_residues=True, apply_ncs=True, init_time='4', nreflections=40000, alphahelix=True,
                 logger=None, resolution=2.0, silent_start=False):

        super(Shelxe, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self._cc_eobs_ecalc = "NA"
        self._cc = "NA"
        self._acl = "NA"
        self._average_cc_delta = "NA"
        self._solution = "NO"
        self._pdbin = pdbin
        self._solvent = solvent
        self._autotracing_ncyc = autotracing_ncyc
        self._density_modif_ncyc = density_modif_ncyc
        self._prune_residues = prune_residues
        self._apply_ncs = apply_ncs
        self._init_time = init_time
        self._nreflections = nreflections
        self._use_f = use_f
        self._resolution = resolution
        self._alphahelix = alphahelix
        self._mtzin = mtzin

    @property
    def wrapper_name(self):
        return "shelxe"

    @property
    def cmd(self):
        return ["shelxe", os.path.basename(self.input_pda)] + [x for x in self.keywords.split()]

    @property
    def keywords(self):
        """Keywords to be used with shelxe"""

        keywrds = '-s%s -a%s -m%s -t%s -l%s' % (self.solvent, self.autotracing_ncyc, self.density_modif_ncyc,
                                                self.init_time, math.ceil(self.nreflections / 10000))
        if self.apply_ncs:
            keywrds += ' -n'
        if self.alphahelix:
            keywrds += ' -q'
        if self.use_f:
            keywrds += ' -f'
        if self.resolution < 2.0:
            keywrds += ' -e1.0'
        if self.prune_residues:
            keywrds += ' -o'

        return keywrds

    @property
    def autotracing_ncyc(self):
        return self._autotracing_ncyc

    @autotracing_ncyc.setter
    def autotracing_ncyc(self, value):
        self._autotracing_ncyc = value

    @property
    def density_modif_ncyc(self):
        return self._density_modif_ncyc

    @density_modif_ncyc.setter
    def density_modif_ncyc(self, value):
        self._density_modif_ncyc = value

    @property
    def prune_residues(self):
        return self._prune_residues

    @prune_residues.setter
    def prune_residues(self, value):
        self._prune_residues = value

    @property
    def apply_ncs(self):
        return self._apply_ncs

    @apply_ncs.setter
    def apply_ncs(self, value):
        self._apply_ncs = value

    @property
    def nreflections(self):
        return self._nreflections

    @nreflections.setter
    def nreflections(self, value):
        self._nreflections = value

    @property
    def alphahelix(self):
        return self._alphahelix

    @alphahelix.setter
    def alphahelix(self, value):
        self._alphahelix = value

    @property
    def use_f(self):
        return self._use_f

    @use_f.setter
    def use_f(self, value):
        self._use_f = value

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, value):
        self._resolution = value

    @property
    def pdbin(self):
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def init_time(self):
        return self._init_time

    @init_time.setter
    def init_time(self, value):
        self._init_time = value

    @property
    def solvent(self):
        return self._solvent

    @solvent.setter
    def solvent(self, value):
        self._solvent = value

    @property
    def cc_eobs_ecalc(self):
        return self._cc_eobs_ecalc

    @cc_eobs_ecalc.setter
    def cc_eobs_ecalc(self, value):
        self._cc_eobs_ecalc = value

    @property
    def cc(self):
        return self._cc

    @cc.setter
    def cc(self, value):
        self._cc = value

    @property
    def acl(self):
        return self._acl

    @acl.setter
    def acl(self, value):
        self._acl = value

    @property
    def mtzin(self):
        return self._mtzin

    @mtzin.setter
    def mtzin(self, value):
        self._mtzin = value

    @property
    def solution(self):
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def average_cc_delta(self):
        return self._average_cc_delta

    @average_cc_delta.setter
    def average_cc_delta(self, value):
        self._average_cc_delta = value

    @property
    def input_pda(self):
        return os.path.join(self.workdir, "shelxe-input.pda")

    @property
    def input_hkl(self):
        return os.path.join(self.workdir, "shelxe-input.hkl")

    @property
    def output_pdb(self):
        return os.path.join(self.workdir, "shelxe-input.pdb")

    @property
    def result_summary(self):
        """String with a summary of the figures of merit obtained"""
        return "Shelxe results: CC - %s   ACL - %s   SOLUTION - %s\n" % (self.cc, self.acl, self.solution)

    def get_scores(self, logfile=None):
        """Extract the figures of merit from the logfile and the pdb output file

        :param logfile: Not in use (default None)
        :type logfile: None
        :returns nothing
        :rtype None
        """

        parser = ShelxeParser(fname=self.output_pdb, stdout=self.logcontents, logger=self.logger)
        parser.parse()
        self.cc, self.acl, self.cc_eobs_ecalc, self.average_cc_delta, self.solution = parser.summary

    def run(self):
        """Run shelxe on the shell"""

        if not self.silent_start:
            self.logger.info(self.wrapper_header)
        self.logger.info('Running %s autotracing cycles with %s cycles of density modification' % (
            self.autotracing_ncyc, self.density_modif_ncyc))

        # Make workdir and get inside
        self.make_workdir()
        tmp_dir = os.getcwd()
        os.chdir(self.workdir)

        # Convert mtzfile into hklfile
        my_mtz2various = Mtz2Various(mtzin=self.mtzin, hklout=self.input_hkl,
                                     workdir=os.path.join(self.workdir, 'mtz2various'))
        my_mtz2various.run()
        my_mtz2various.make_logfile()

        # Run shelxe
        copyfile(self.pdbin, self.input_pda)
        self.logger.debug(" ".join(self.cmd))
        self.logcontents = cexec(self.cmd, stdin=self.keywords, permit_nonzero=True)

        # Get the scores
        self.get_scores()
        self.logger.info(self.result_summary)

        os.chdir(tmp_dir)
