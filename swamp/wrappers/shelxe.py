import os
import math
from pyjob import cexec
from shutil import copyfile
from swamp.parsers import ShelxeParser
from swamp.wrappers.wrapper import Wrapper
from swamp.wrappers.mtz2various import Mtz2Various


class Shelxe(Wrapper):
    """Wrapper around shelxe

    :param str workdir: working directory
    :param str mtzin: target structure mtz file name
    :param str pdbin: pdb file name with the refined search model
    :param str solvent: estimated solvent content for the crystal (default 0.45)
    :param str autotracing_ncyc: number of autotracing cycles to run
    :param str density_modif_ncyc: number of density modification cycles to run
    :param bool bool use_f: if True use F labels in the mtz file
    :param bool prune_residues: if True set -o option (default True)
    :param bool apply_ncs: sets parameter -n (default True)
    :param str init_time: parameter to set -t (default True)
    :param int nreflections: number of reflections in the mtz file (default 40000)
    :param bool alphahelix: if True sets the -t parameter (default True)
    :param float resolution: resolution of the data in the mtz file (default 2.0)
    :param bool silent_start: if True, the logger will not display the start banner (default False)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :ivar bool error: if True an error has occurred along the process
    :ivar str cc_eobs_ecalc: the correlation coeff. between Eobs and Ecalc
    :ivar str cc: correlation coeff. obtained with the best tracing cycle
    :ivar str acl: average chain length obtained with the best tracing cycle
    :ivar str average_cc_delta: the average delta of correlation coeff. between tracing cycles
    :ivar str solution: 'YES' if correlation coeff > 25 otherwise 'NO'

    :examples:

    >>> from swamp.wrappers import Shelxe
    >>> my_shelxe = Shelxe('<workdir>', '<mtzin>', '<pdbin>')
    >>> my_shelxe.run()
    >>> my_shelxe.make_logfile()

    """

    def __init__(self, workdir, pdbin, mtzin, solvent='0.45', autotracing_ncyc='20', density_modif_ncyc='10',
                 use_f=True, prune_residues=True, apply_ncs=True, init_time='4', nreflections=40000, alphahelix=True,
                 logger=None, resolution=2.0, silent_start=False):

        super(Shelxe, self).__init__(workdir=workdir, logger=logger, silent_start=silent_start)

        self.cc_eobs_ecalc = "NA"
        self.cc = "NA"
        self.acl = "NA"
        self.average_cc_delta = "NA"
        self.solution = "NO"
        self.pdbin = pdbin
        self.solvent = solvent
        self.autotracing_ncyc = autotracing_ncyc
        self.density_modif_ncyc = density_modif_ncyc
        self.prune_residues = prune_residues
        self.apply_ncs = apply_ncs
        self.init_time = init_time
        self.nreflections = nreflections
        self.use_f = use_f
        self.resolution = resolution
        self.alphahelix = alphahelix
        self.mtzin = mtzin

    @property
    def wrapper_name(self):
        """The name of this wrapper (shelxe)"""
        return "shelxe"

    @property
    def cmd(self):
        """Command to be executed on the shell"""
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
    def input_pda(self):
        """The input .pda file name for shelxe"""
        return os.path.join(self.workdir, "shelxe-input.pda")

    @property
    def input_hkl(self):
        """The input .hkl file name for shelxe"""
        return os.path.join(self.workdir, "shelxe-input.hkl")

    @property
    def output_pdb(self):
        """The output .pdb file name created by shelxe"""
        return os.path.join(self.workdir, "shelxe-input.pdb")

    @property
    def summary_results(self):
        """String representation of a summary of the figures of merit obtained"""
        return "Shelxe results: CC - %s   ACL - %s   SOLUTION - %s\n" % (self.cc, self.acl, self.solution)

    def get_scores(self, logfile=None):
        """Extract the figures of merit from the logfile and the pdb output file using \
        :py:obj:`~swamp.parsers.shelxeparser.ShelxeParser`

        :param None logfile: Not in use (default None)
        """

        parser = ShelxeParser(fname=self.output_pdb, stdout=self.logcontents, logger=self.logger)
        parser.parse()
        self.cc, self.acl, self.cc_eobs_ecalc, self.average_cc_delta, self.solution = parser.summary

    def run(self):
        """Run :py:attr:`~swamp.wrappers.shelxe.Shelxe.cmd` and store the results"""

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
