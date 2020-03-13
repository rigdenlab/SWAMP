import os
import shutil
import collections
from swamp.wrappers.wrapper import Wrapper
from swamp.wrappers.wrefmac import wRefmac
from swamp.parsers import MtzParser, PhaserParser
from phaser import InputMR_DAT, runMR_DAT, InputMR_AUTO, runMR_AUTO


class Phaser(Wrapper):
    """Phaser wrapper

    Class with methods to execute phaser and datas tructures to contain the obtained results

    :param str workdir: working directory
    :param srt mtzfile: target structure mtz file name
    :param int mw: molecular weigth of the target structure
    :param int nchains_asu: number of chains found in the assymetric unit (default 1)
    :param str sgalternative: parameter to be passed to phaser as \
    :py:func:`phaser.InputMR_AUTO.SGAL_SELE` (default 'NONE')
    :param bool early_kill: if True will send signal to abort the run if the figures of merit are to low (default True)
    :param int timeout: set a limit kill time for phaser execution (default 360)
    :param int threads: number of threads to be used with phaser (default 1)
    :param packcutoff: parameter to be passed to phaser as :py:func:`phaser.InputMR_AUTO.setPACK_CUTO` (default None)
    :type packcutoff: float, None
    :param peaks_rotcutoff: parameter to be passed to phaser as :py:func:`phaser.InputMR_AUTO.setPEAK_ROTA_CUTO`\
     (default None)
    :type peaks_rotcutoff: float, None
    :param str phased_mtz: target's mtz filename containing phases (default: None)
    :param bool silent_start: if True, the logger will not display the start banner (default False)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the wrapper (default None)
    :ivar str LLG: LLG of the placed search model
    :ivar str TFZ: TFZ of the placed search model
    :ivar str RFZ: RFZ of the placed search model
    :ivar str VRMS: VRMS calculated for the placed search model
    :ivar str eLLG: eLLG calculated for the search model
    :ivar bool error: if True an error has occurred along the process
    :ivar bool abort_suggested: if True the figures of merit are low and it is suggested to abort \
    :py:obj:`~swamp.mr.mrrun.MrRun`
    :ivar str ovarall_CC: overall correlation coefficient between phases at \
    :py:attr:`~swamp.wrappers.wphaser.Phaser.phased_mtz` and the placed search model
    :ivar str local_CC: local correlation coefficient between phases at \
    :py:attr:`~swamp.wrappers.wphaser.Phaser.phased_mtz` and the placed search model
    :ivar tuple solution: a named tuple with information about an already existing solution
    :ivar list searchmodel_ids_list: a list with the search model identifiers
    :ivar dict searchmodel_dict: a dictionary with the search model information

    :example:

    >>> from swamp.wrappers import Phaser
    >>> my_phaser = Phaser('<workdir>', '<mtzfile>', <mw>)
    >>> my_phaser.add_searchmodel('<pdbfile>')
    >>> my_phaser.run()
    >>> my_phaser.make_logfile()
    """

    def __init__(self, workdir, mtzfile, mw, nchains_asu=1, sgalternative='NONE', logger=None, early_kill=True,
                 timeout=360, threads=1, phased_mtz=None, silent_start=False, packcutoff=None, peaks_rotcutoff=None):

        super(Phaser, self).__init__(logger=logger, workdir=workdir, silent_start=silent_start)

        self.abort_suggested = False
        self.RFZ = "NA"
        self.TFZ = "NA"
        self.LLG = "NA"
        self.eLLG = "NA"
        self.VRMS = "NA"
        self.local_CC = "NA"
        self.overall_CC = "NA"
        self.phased_mtz = phased_mtz
        self.input_mr_dat = None
        self.input_mr_auto = None
        self.run_mr_auto = None
        self.run_mr_dat = None
        self.early_kill = early_kill
        self.threads = threads
        self.timeout = timeout
        self.solution = None
        self.mtzfile = mtzfile
        self.sgalternative = sgalternative
        self.nchains_asu = nchains_asu
        self.mw = mw
        self.packcutoff = packcutoff
        self.peaks_rotcutoff = peaks_rotcutoff
        self.searchmodel_dict = {}
        self.searchmodel_ids_list = []

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['input_mr_dat']
        del d['input_mr_auto']
        del d['run_mr_auto']
        del d['run_mr_dat']
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)

    # ------------------ Properties ------------------

    @property
    def keywords(self):
        """No keywords passed to phaser through stdin"""
        return None

    @property
    def cmd(self):
        """No command run through the shell to execute phaser"""
        return None

    @property
    def wrapper_name(self):
        """The name of this `~swamp.wrapper.wrapper.Wrapper` child class (phaser)"""
        return "phaser"

    @property
    def logfile(self):
        """Log file name with the logging contents of phaser"""
        return os.path.join(self.workdir, 'phaser_PDB_{}_out.log')

    @property
    def summary_results(self):
        """A string representation of a summary of the figures of merit"""
        return "Phaser results: LLG - %s   TFZ - %s   Local CC - %s   Overall CC - %s" % (
            self.LLG, self.TFZ, self.local_CC, self.overall_CC)

    @property
    def pdbout_tmp(self):
        """A temporary pdb output file name template for placed search models"""
        return os.path.join(self.workdir, 'phaser_PDB_{}_out.pdb')

    @property
    def refined_solutions_dir(self):
        """A directory where placed search models will be refined"""
        return os.path.join(self.workdir, 'refined_solutions')

    @property
    def refmac_pdbout(self):
        """A pdb output file name template for refined search models"""
        return os.path.join(self.refined_solutions_dir, 'refined_PDB_{}_out.pdb')

    @property
    def refmac_logfile(self):
        """A log file name template for placed search models"""
        return os.path.join(self.refined_solutions_dir, 'refined_PDB_{}_out.log')

    @property
    def root(self):
        """Root file name to be used with phaser"""
        return "%s_phaser" % os.path.basename(self.mtzfile[:-4])

    @property
    def top_PDBfile(self):
        """PDB file name corresponding with the top solution"""
        return "%s.1.pdb" % os.path.basename(self.root)

    @property
    def top_MTZfile(self):
        """MTZ file name corresponding with the top solution"""
        return "%s.1.mtz" % os.path.basename(self.root)

    @property
    def top_searchmodel(self):
        """First search model at :py:attr:`~swamp.wrappers.wphaser.Phaser.searchmodel_ids_list`"""
        if len(self.searchmodel_ids_list) == 0:
            return None
        else:
            return self.searchmodel_dict[self.searchmodel_ids_list[0]]

    @property
    def solution_infotemplate(self):
        """Named tuple to serve as a template to contain information about the fixed solution (if any is provided)"""
        return collections.namedtuple("solution", ["sol_fname", "pdbfile", "ermsd", "ident"])

    @property
    def searchmodel_infotemplate(self):
        """Named tuple to serve as a template to contain information about the search models"""
        return collections.namedtuple("SearchModelInfo", ["pdbfile", "ermsd", "nsearch", "disable_check", "id"])

    # ------------------ Methods ------------------

    def _run(self):
        """Override :py:func:`~swamp.wrappers.wrapper.Wrapper._run`"""
        pass

    def register_solution(self, pdbfile, ermsd=None, ident=None, sol_fname=None):
        """Add information for an already existing solution

        :param str pdbfile: the pdb file name with the existing solution
        :param float ermsd: the eRMSD used to obtain this particular solution
        :param float ident: the identity used to obtain this particular solution
        :param str sol_fname: the .sol file for this solution
        """

        self.solution = self.solution_infotemplate(pdbfile=pdbfile, ermsd=ermsd, ident=ident, sol_fname=sol_fname)

    def add_searchmodel(self, pdbfile, id, ermsd=0.1, nsearch=1, disable_check=True):
        """Add a search model to the phaser run

        :param str pdbfile: the pdb file name with the search model
        :param str id: unique identifier for this search model
        :param float ermsd: the eRMSD to be used with this search model
        :param int nsearch: number of copies to be searched
        :param bool disable_check: argument passed to :py:obj:`phaser.InputMR_AUTO.setENSE_DISA_CHEC()`
        """

        if not os.path.isfile(pdbfile):
            self.logger.error('Searchmodel pdb file was not found! %s' % pdbfile)
            self.error = True
            return

        if id in self.searchmodel_dict.keys():
            self.logger.error('Searchmodel with the same ID (%s) has already been added!' % id)
            self.error = True
            return

        self.searchmodel_ids_list.append(id)
        self.searchmodel_dict[id] = self.searchmodel_infotemplate(pdbfile=pdbfile, ermsd=ermsd, nsearch=nsearch, id=id,
                                                                  disable_check=disable_check)

    def load_mr_dat(self, mute=True):
        """Load the input data for the MR run into :py:obj:`phaser.InputMR_DAT()`

        :param bool mute: argument passed to :py:obj:`phaser.InputMR_DAT.setMUTE()`
        """

        self.input_mr_dat = InputMR_DAT()
        self.input_mr_dat.setHKLI(self.mtzfile)

        # Set reflection data fields (I and SIGI preferred over FP and SIGFP)
        mtz_head = MtzParser(self.mtzfile)
        mtz_head.parse()
        if mtz_head.i is not None and mtz_head.sigi is not None:
            self.input_mr_dat.setLABI_I_SIGI(mtz_head.i, mtz_head.sigi)
        elif mtz_head.f is not None and mtz_head.sigf is not None:
            self.input_mr_dat.setLABI_F_SIGF(mtz_head.f, mtz_head.sigf)
        elif mtz_head.iplus is not None and mtz_head.sigiplus is not None:
            self.input_mr_dat.setLABI_I_SIGI(mtz_head.iplus, mtz_head.sigiplus)
        elif mtz_head.fplus is not None and mtz_head.sigfplus is not None:
            self.input_mr_dat.setLABI_F_SIGF(mtz_head.fplus, mtz_head.sigfplus)
        elif mtz_head.iminus is not None and mtz_head.sigiminus is not None:
            self.input_mr_dat.setLABI_I_SIGI(mtz_head.iminus, mtz_head.sigiminus)
        elif mtz_head.fminus is not None and mtz_head.sigfminus is not None:
            self.input_mr_dat.setLABI_F_SIGF(mtz_head.fminus, mtz_head.sigfminus)
        else:
            self.logger.error("Cannot find mtz columns to use as input for phaser!")
            self.error = True
            return

        self.input_mr_dat.setSGAL_SELE(self.sgalternative)
        self.input_mr_dat.setMUTE(mute)

        # Load the data into the mr run
        self.run_mr_dat = runMR_DAT(self.input_mr_dat)
        if not self.run_mr_dat.Success():
            self.logger.error("Phaser didn't register success!")
            self.error = True

    def run_auto_mr(self, searchmodel_id, mute=True, xyzout=True, keywrds=True):
        """Make a method to run phaser AutoMR mode using :object:`Phaser.runMR_AUTO()`

        :param str searchmodel_id: key in the :py:attr:`~swamp.wrappers.wphaser.Phaser.searchmodel_dict` that \
        corresponds with the search model to be placed
        :param bool mute: argument passed to :py:obj:`phaser.InputMR_AUTO.setMUTE()`
        :param bool xyzout: argument passed to :py:obj:`phaser.InputMR_AUTO.setXYZO()`
        :param bool keywrds: argument passed to :py:obj:`phaser.InputMR_AUTO.setKEYW()`
        """

        if searchmodel_id not in self.searchmodel_dict.keys():
            self.logger.error('Unable to find requested searchmodel with ID %s!' % searchmodel_id)
            self.error = True
            return

        self.logger.info('Running auto MR for search model %s' % searchmodel_id)

        self.input_mr_auto = InputMR_AUTO()
        searchmodel = self.searchmodel_dict[searchmodel_id]

        # General parameters
        self.input_mr_auto.setJOBS(self.threads)
        self.input_mr_auto.setHKLI(self.mtzfile)
        self.input_mr_auto.setXYZO(xyzout)
        self.input_mr_auto.setROOT(self.root)
        self.input_mr_auto.setKEYW(keywrds)
        self.input_mr_auto.setKILL_TIME(self.timeout)
        self.input_mr_auto.setREFL_DATA(self.run_mr_dat.getREFL_DATA())
        self.input_mr_auto.setMUTE(mute)
        self.input_mr_auto.setSGAL_SELE(self.sgalternative)
        if self.packcutoff is not None:
            self.input_mr_auto.setPACK_SELE('PERCENT')
            self.input_mr_auto.setPACK_CUTO(self.packcutoff)
        if self.peaks_rotcutoff is not None:
            self.input_mr_auto.setPEAK_ROTA_SELE('PERCENT')
            self.input_mr_auto.setPEAK_ROTA_CUTO(self.peaks_rotcutoff)
        self.input_mr_auto.setCOMP_BY("ASU")
        self.input_mr_auto.addCOMP_PROT_MW_NUM(self.mw, self.nchains_asu)

        # Load the search models
        self.input_mr_auto.addENSE_PDB_RMS("PDB_%s" % searchmodel_id, searchmodel.pdbfile, float(searchmodel.ermsd))
        self.input_mr_auto.setENSE_DISA_CHEC("PDB_%s" % searchmodel_id, searchmodel.disable_check)
        self.input_mr_auto.addSEAR_ENSE_NUM("PDB_%s" % searchmodel_id, searchmodel.nsearch)

        # Load fixed solution if any
        if self.solution is not None:
            if self.solution.ermsd is not None:
                self.input_mr_auto.addENSE_PDB_RMS("FIX", self.solution.pdbfile, float(self.solution.ermsd))
            elif self.solution.ident is not None:
                self.input_mr_auto.addENSE_PDB_ID("FIX", self.solution.pdbfile, float(self.solution.ident))
            else:
                self.logger.error("NEED TO PROVIDE AT LEAST RMSD OR IDENT FOR THE FIXED SOLUTION!")
                self.error = True
                return
            self.input_mr_auto.addSOLU_ORIG_ENSE("FIX")

        # Run mr AUTO
        self.run_mr_auto = runMR_AUTO(self.input_mr_auto)
        self.logcontents = self.run_mr_auto.summary()
        # Check for errors
        self._check_error()
        if not self.error:
            shutil.move(self.top_PDBfile, self.pdbout_tmp.format(searchmodel_id))
            shutil.copyfile(self.pdbout_tmp.format(searchmodel_id), self.pdbout)
            self.get_scores(self.pdbout_tmp.format(searchmodel_id))
            self.logger.info(self.summary_results)

    def run(self):
        """Run phaser

        Load the MR data and run mode MR_AUTO with each of the provided search models, intercalating refinement cycles
        """

        if not self.silent_start:
            self.logger.info(self.wrapper_header)

        if len(self.searchmodel_dict.keys()) == 0:
            self.logger.error('No search model found!')
            return

        self.make_workdir()
        tmp_dir = os.getcwd()
        os.chdir(self.workdir)

        # Load the MR data
        self.logger.debug("Phaser: loading input for MR_DAT")
        self.load_mr_dat()
        if self.error:
            return

        # Run MR auto incrementally with provided search models
        self.logger.info("Running incremental MR_AUTO with %s searchmodels" % len(self.searchmodel_ids_list))
        for searchmodel_id in self.searchmodel_ids_list:
            self.run_auto_mr(searchmodel_id)
            self.make_logfile(fname=self.logfile.format(searchmodel_id))
            if self.error:
                self.logger.warning('Previous error prevents phaser to proceed with the search')
                break
            if searchmodel_id == self.searchmodel_ids_list[-1]:
                self.logger.info('Finished incremental search\n')
                break
            self.refine_solution(searchmodel_id)
            if self.error:
                self.logger.warning('Previous error during solution refinement prevents phaser to proceed'
                                    ' with the search')
                break
            self.register_solution(pdbfile=self.refmac_pdbout.format(searchmodel_id),
                                   ermsd=self.searchmodel_dict[searchmodel_id].ermsd)

        os.chdir(tmp_dir)

    def get_scores(self, logfile=None):
        """Get the figures of merit obtained with phaser using a :py:obj:`~swamp.parsers.phaserparser.PhaserParser` \
        instance

        :param logfile: log's file name where the figures of merit can be found (default None)
        :type logfile: None, str
        """

        parser = PhaserParser(fname=logfile, stdout=self.logcontents, logger=self.logger)
        parser.parse()
        self.LLG, self.TFZ, self.RFZ, self.eLLG, self.VRMS = parser.summary

        if parser.error:
            self.error = True
            self.logger.warning('Previous errors while parsing phaser output detected!')
            return

        if self.phased_mtz is not None:
            self.local_CC, self.overall_CC = self.get_cc(os.path.join(self.workdir, "phenix"), self.phased_mtz, logfile)
            if self.early_kill and (self.local_CC == "NA" or float(self.local_CC) < 0.2):
                self.abort_suggested = True
            else:
                self.abort_suggested = False

    def refine_solution(self, searchmodel_id):
        """Refine the placed search model using :py:obj:`~swamp.wrappers.wrefmac.wRefmac`

        :param str searchmodel_id: key in the :py:attr:`~swamp.wrappers.wphaser.Phaser.searchmodel_dict` that \
        corresponds with the search model to be refined
        """

        self.logger.info('Refining solution')

        refmac = wRefmac(workdir=self.refined_solutions_dir, mtzin=self.mtzfile, logger=self.logger, silent_start=True,
                         pdbin=self.pdbout_tmp.format(searchmodel_id), phased_mtz=self.phased_mtz)
        refmac.run()
        refmac.make_logfile(fname=self.refmac_logfile.format(searchmodel_id))

        # Check for errors and copy the file to the corresponding pdb output
        if refmac.error or not os.path.isfile(refmac.pdbout):
            self.error = True
            self.logger.warning('An error occurred during refinement!')
        else:
            shutil.copyfile(refmac.pdbout, self.refmac_pdbout.format(searchmodel_id))

    # ------------------ Static and hidden methods ------------------

    def _check_error(self):
        """Check for errors in the :object:`Phaser.runMR_AUTO()` object

        Sets the :py:attr:`~swamp.wrappers.wphaser.Phaser.error` to True is an error is found after \
        :py:func:`~swamp.wrappers.wphaser.Phaser.run_auto_mr`
        """

        if not self.run_mr_auto.Success():
            if self.run_mr_auto.ErrorName() == 'KILL-TIME ELAPSED':
                self.logger.warning("Phaser exceeded the kill-time! (%s min)" % self.timeout)
            else:
                self.logger.warning("Phaser did not register success!")
                self.error = True
                return

        if not self.run_mr_auto.foundSolutions():
            self.logger.warning("Phaser did not find a solution!")
            self.error = True
            return

        if not os.path.isfile(self.top_PDBfile):
            self.logger.warning("Phaser did not produce a pdb output file!")
            self.error = True
            return

        if not os.path.isfile(self.top_MTZfile):
            self.logger.warning("Phaser did not produce a mtz output file!")
            return
