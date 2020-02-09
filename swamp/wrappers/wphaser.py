import os
import shutil
import collections
from simbad.util import mtz_util
from swamp.wrappers.wrapper import Wrapper
from swamp.wrappers.wrefmac import wRefmac
from swamp.mr.searchmodel_prepare.polyala import PolyALA
from phaser import InputMR_DAT, runMR_DAT, InputMR_AUTO, runMR_AUTO


class Phaser(Wrapper):
    """Phaser wrapper

    Class with methods to execute phaser and datas tructures to contain the obtained results

     :param workdir: working directory
     :type workdir: str, None
     :param mtzfile: target structure mtz file name
     :type mtzfile: str
     :param mw: molecular weigth of the target structure
     :type mw: int
     :param nchains_asu: number of chains found in the assymetric unit (default 1)
     :type nchains_asu: int
     :param sgalternative: parameter to be passed to phaser as SGAL_SELE (default 'NONE')
     :type sgalternative: str
     :param early_kill: if True will send signal to abort the run if the figures of merit are to low (default True)
     :type early_kill: bool
     :param timeout: set a limit kill time for phaser execution (default 360)
     :type timeout: int
     :param threads: number of threads to be used with phasern(default 1)
     :type threads: int
     :param packcutoff: parameter to be passed to phaser as setPACK_CUTO (default None)
     :type packcutoff: float, None
     :param peaks_rotcutoff: parameter to be passed to phaser as setPEAK_ROTA_CUTO (default None)
     :type peaks_rotcutoff: float, None
     :param phased_mtz: target's mtz filename containing phases (default: None)
     :type phased_mtz: str
     :param silent_start: if True, the logger will not display the start banner (default False)
     :type silent_start: bool
     :param logger: logging interface to be used (default None)
     :type logger: None, :object:`swamp.logger.swamplogger.SwampLogger`
     :ivar LLG: LLG of the placed search model
     :type LLG: str
     :ivar TFZ: TFZ of the placed search model
     :type TFZ: str
     :ivar RFZ: RFZ of the placed search model
     :type RFZ: str
     :ivar VRMS: VRMS calculated for the placed search model
     :type VRMS: str
     :ivar eLLG: eLLG calculated for the search model
     :type eLLG: str
     :ivar error: if True an error has occurred along the process
     :type error: bool

     :examples

    >>> from swamp.wrappers.wphaser import Phaser
    >>> my_phaser = Phaser('<workdir>', '<mtzfile>', <mw>)
    >>> my_phaser.add_searchmodel('<pdbfile>')
    >>> my_phaser.run()
    >>> my_phaser.make_logfile()
    """

    def __init__(self, workdir, mtzfile, mw, nchains_asu=1, sgalternative='NONE', logger=None, early_kill=True,
                 timeout=360, threads=1, phased_mtz=None, silent_start=False, packcutoff=None, peaks_rotcutoff=None):

        super(Phaser, self).__init__(logger=logger, workdir=workdir, silent_start=silent_start)

        self._failure = None
        self._abort_suggested = False
        self._RFZ = "NA"
        self._TFZ = "NA"
        self._LLG = "NA"
        self._eLLG = "NA"
        self._VRMS = "NA"
        self._local_CC = "NA"
        self._overall_CC = "NA"
        self._phased_mtz = phased_mtz
        self._input_mr_dat = None
        self._input_mr_auto = None
        self._run_mr_auto = None
        self._run_mr_dat = None
        self._early_kill = early_kill
        self._threads = threads
        self._timeout = timeout
        self._solution = None
        self._mtzfile = mtzfile
        self._sgalternative = sgalternative
        self._nchains_asu = nchains_asu
        self._mw = mw
        self._packcutoff = packcutoff
        self._peaks_rotcutoff = peaks_rotcutoff
        self._searchmodel_dict = {}
        self._searchmodel_ids_list = []

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['_input_mr_dat']
        del d['_input_mr_auto']
        del d['_run_mr_auto']
        del d['_run_mr_dat']
        return d

    def __setstate__(self, d):
        self.__dict__.update(d)

    # ------------------ Properties ------------------

    @property
    def keywords(self):
        return None

    @property
    def cmd(self):
        return None

    @property
    def wrapper_name(self):
        return "phaser"

    @property
    def logfile(self):
        return os.path.join(self.workdir, 'phaser_PDB_{}_out.log')

    @property
    def summary_results(self):
        return "Phaser results: LLG - %s   TFZ - %s   Local CC - %s   Overall CC - %s" % (
            self.LLG, self.TFZ, self.local_CC, self.overall_CC)

    @property
    def pdbout_tmp(self):
        return os.path.join(self.workdir, 'phaser_PDB_{}_out.pdb')

    @property
    def mtzfile(self):
        return self._mtzfile

    @mtzfile.setter
    def mtzfile(self, value):
        self._mtzfile = value

    @property
    def packcutoff(self):
        return self._packcutoff

    @packcutoff.setter
    def packcutoff(self, value):
        self._packcutoff = value

    @property
    def peaks_rotcutoff(self):
        return self._peaks_rotcutoff

    @peaks_rotcutoff.setter
    def peaks_rotcutoff(self, value):
        self._peaks_rotcutoff = value

    @property
    def mw(self):
        return self._mw

    @mw.setter
    def mw(self, value):
        self._mw = value

    @property
    def sgalternative(self):
        return self._sgalternative

    @sgalternative.setter
    def sgalternative(self, value):
        self._sgalternative = value

    @property
    def nchains_asu(self):
        return self._nchains_asu

    @nchains_asu.setter
    def nchains_asu(self, value):
        self._nchains_asu = value

    @property
    def refined_solutions_dir(self):
        return os.path.join(self.workdir, 'refined_solutions')

    @property
    def refmac_pdbout(self):
        return os.path.join(self.refined_solutions_dir, 'refined_PDB_{}_out.pdb')

    @property
    def refmac_logfile(self):
        return os.path.join(self.refined_solutions_dir, 'refined_PDB_{}_out.log')

    @property
    def phased_mtz(self):
        return self._phased_mtz

    @phased_mtz.setter
    def phased_mtz(self, value):
        self._phased_mtz = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        self._threads = value

    @property
    def timeout(self):
        return self._timeout

    @timeout.setter
    def timeout(self, value):
        self._timeout = value

    @property
    def early_kill(self):
        return self._early_kill

    @early_kill.setter
    def early_kill(self, value):
        self._early_kill = value

    @property
    def input_mr_dat(self):
        return self._input_mr_dat

    @input_mr_dat.setter
    def input_mr_dat(self, value):
        self._input_mr_dat = value

    @property
    def input_mr_auto(self):
        return self._input_mr_auto

    @input_mr_auto.setter
    def input_mr_auto(self, value):
        self._input_mr_auto = value

    @property
    def run_mr_auto(self):
        return self._run_mr_auto

    @run_mr_auto.setter
    def run_mr_auto(self, value):
        self._run_mr_auto = value

    @property
    def run_mr_dat(self):
        return self._run_mr_dat

    @run_mr_dat.setter
    def run_mr_dat(self, value):
        self._run_mr_dat = value

    @property
    def solution(self):
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def root(self):
        return "%s_phaser" % os.path.basename(self.mtzfile[:-4])

    @property
    def abort_suggested(self):
        """True if the figures of merit are low and the search model is likely to be incorrectly placed"""
        return self._abort_suggested

    @abort_suggested.setter
    def abort_suggested(self, value):
        self._abort_suggested = value

    @property
    def RFZ(self):
        return self._RFZ

    @RFZ.setter
    def RFZ(self, value):
        self._RFZ = value

    @property
    def TFZ(self):
        return self._TFZ

    @TFZ.setter
    def TFZ(self, value):
        self._TFZ = value

    @property
    def LLG(self):
        return self._LLG

    @LLG.setter
    def LLG(self, value):
        self._LLG = value

    @property
    def eLLG(self):
        return self._eLLG

    @eLLG.setter
    def eLLG(self, value):
        self._eLLG = value

    @property
    def VRMS(self):
        return self._VRMS

    @VRMS.setter
    def VRMS(self, value):
        self._VRMS = value

    @property
    def local_CC(self):
        """local CC obtained with this phaser run (not None only if pdb benchmark file is provided)"""
        return self._local_CC

    @local_CC.setter
    def local_CC(self, value):
        self._local_CC = value

    @property
    def overall_CC(self):
        """overall CC obtained with this phaser run (not None only if pdb benchmark file is provided)"""
        return self._overall_CC

    @overall_CC.setter
    def overall_CC(self, value):
        self._overall_CC = value

    @property
    def searchmodel_dict(self):
        """Dictionary of search models to be used in the run, keys correspond with their unique id"""
        return self._searchmodel_dict

    @searchmodel_dict.setter
    def searchmodel_dict(self, value):
        self._searchmodel_dict = value

    @property
    def searchmodel_ids_list(self):
        """List of search models to be used in the phaser run"""
        return self._searchmodel_ids_list

    @searchmodel_ids_list.setter
    def searchmodel_ids_list(self, value):
        self._searchmodel_ids_list = value

    @property
    def top_searchmodel(self):
        """First search model of the list"""
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

    def register_solution(self, pdbfile, ermsd=None, ident=None, sol_fname=None):
        """Add information for an already existing solution

        :param pdbfile: the pdb file name with the existing solution
        :type pdbfile: str
        :param ermsd: the eRMSD used to obtain this particular solution
        :type ermsd: float
        :param ident: the identity used to obtain this particular solution
        :type ident: float
        :param sol_fname: the .sol file for this solution
        :type sol_fname: str
        :returns nothing
        :rtype None
        """

        self.solution = self.solution_infotemplate(pdbfile=pdbfile, ermsd=ermsd, ident=ident, sol_fname=sol_fname)

    def add_searchmodel(self, pdbfile, id, ermsd=0.1, nsearch=1, disable_check=True):
        """Add a search model to the phaser run

        :param pdbfile: the pdb file name with the search model
        :type pdbfile: str
        :param id: unique identifier for this search model
        :type id: str, int
        :param ermsd: the eRMSD to be used with this search model
        :type ermsd: float
        :param nsearch: number of copies to be searched
        :type nsearch: int
        :param disable_check: argument passed to :object:`Phaser.InputMR_AUTO.setENSE_DISA_CHEC()`
        :type disable_check: bool
        :returns nothing
        :rtype None
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
        """Load the input data for the MR run into :object:`Phaser.InputMR_DAT()`

        :param mute: argument passed to :object:`Phaser.InputMR_DAT.setMUTE()`
        :type mute: bool
        :returns nothing
        :rtype None
        """

        self.input_mr_dat = InputMR_DAT()
        self.input_mr_dat.setHKLI(self.mtzfile)

        # Set reflection data fields (I and SIGI preferred over FP and SIGFP)
        mtz_head = mtz_util.GetLabels(self.mtzfile)
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

        :param searchmodel_id: key in the searchmodel_dict attribute that corresponds with the search model to be placed
        :type searchmodel_id: str, int
        :param mute: argument passed to :object:`Phaser.InputMR_AUTO.setMUTE()`
        :type mute: bool
        :param xyzout: argument passed to :object:`Phaser.InputMR_AUTO.setXYZO()`
        :type xyzout: bool
        :param keywrds: argument passed to :object:`Phaser.InputMR_AUTO.setKEYW()`
        :type keywrds: bool
        :returns nothing
        :rtype None
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
            shutil.move(self.run_mr_auto.getTopPdbFile(), self.pdbout_tmp.format(searchmodel_id))
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
        """ Get the figures of merit obtained with phaser

        :param logfile: log's file name where the figures of merit can be found (default None)
        :type logfile: None, str
        :returns nothing
        :rtype None
        """

        self.LLG, self.TFZ, self.RFZ = self._parse_pdbout(pdbout=logfile)
        self.eLLG, self.VRMS = self._parse_logfile(self.logcontents)
        if self.LLG == "NA" or self.TFZ == "NA":
            self.logger.error("Unable to parse TFZ (%s) and LLG (%s)" % (self.TFZ, self.LLG))
            self.error = True
            return
        if self.phased_mtz is not None:
            self.local_CC, self.overall_CC = self.get_cc(os.path.join(self.workdir, "phenix"), self.phased_mtz, logfile)
            if self.early_kill and (self.local_CC == "NA" or float(self.local_CC) < 0.1):
                self.abort_suggested = True
            else:
                self.abort_suggested = False

    def refine_solution(self, searchmodel_id):
        """Refine the solution for a given searchmodel using refmac5

        :param searchmodel_id: key in the searchmodel_dict attribute that corresponds with the search model to be placed
        :type searchmodel_id: str, int
        :returns nothing
        :rtype None
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

        Sets the error attribute to True is an error is found after running MR auto

        :returns nothing
        :rtype None
        """

        if not self.run_mr_auto.Success():
            self.logger.error("Phaser did not register success!")
            self.error = True

        elif not os.path.isfile(self.run_mr_auto.getTopPdbFile()):
            self.logger.error("Phaser did not produce a pdb output file!")
            self.error = True

        elif not self.run_mr_auto.foundSolutions():
            self.logger.error("Phaser did not find a solution!")
            self.error = True

        elif not os.path.isfile(self.run_mr_auto.getTopMtzFile()):
            self.logger.warning("Phaser did not produce a mtz file!")

        elif self.run_mr_auto.Failure():
            self.logger.error("Phaser registered failure!")
            self.error = True

    @staticmethod
    def remove_side_chains(pdbfile):
        """Method to remove the side chains of a given pdb file

        :param pdbfile: filename of the pdb file to be modified
        :type pdbfile: str
        :returns no value
        :rtype None
        """

        tmp_file = os.path.join(os.path.dirname(pdbfile), "phaser_out_polyALA.pdb")
        PolyALA.truncate_polyALA(pdbin=pdbfile, pdbout=tmp_file)
        PolyALA.renumber_residues(pdbin=tmp_file)
        PolyALA.transfer_flags_pdb(pdb_ref=pdbfile, pdb_file=tmp_file)
        shutil.move(tmp_file, pdbfile)

    @staticmethod
    def check_intensities(input_mtz):
        """Determine if a mtz file contains intensities (preferred over amplitudes)

        :param input_mtz: mtz file name of interest
        :type input_mtz: str
        :returns True if intensity labels are found in the mtz file
        :rtype bool
        """

        mtz_head = mtz_util.GetLabels(input_mtz)
        if mtz_head.i is not None and mtz_head.sigi is not None:
            return True
        else:
            return False

    @staticmethod
    def _parse_pdbout(pdbout):
        """Method to get LLG, TFZ and RFZ out of a pdbout

        :param pdbout: the file name of the pdb output generated by phaser
        :type pdbout: str
        :returns a tuple with the LLG (str), the TFZ (str) and the RFZ (str)
        :rtype: tuple
        """

        LLG = "NA"
        TFZ = "NA"
        RFZ = "NA"
        # Parse the pdbout for LLG, TFZ and RFZ
        with open(pdbout, "r") as fhandle:
            for line in fhandle:
                if "Log-Likelihood Gain" in line:
                    LLG = line.split(":")[1].rstrip().lstrip()
                    continue
                elif "TFZ" in line or "RFZ" in line:
                    for score in line.split():
                        if "TFZ==" in score:
                            TFZ = score.split("==")[1].rstrip().lstrip()
                        elif "TFZ=" in score:
                            TFZ = score.split("=")[1].rstrip().lstrip()
                        elif "RFZ==" in score:
                            RFZ = score.split("==")[1].rstrip().lstrip()
                        elif "RFZ=" in score:
                            RFZ = score.split("=")[1].rstrip().lstrip()
                    continue
                # If LLG and TFZ are stored, break
                if TFZ != "NA" and LLG != "NA" and RFZ != "NA":
                    break
        return LLG, TFZ, RFZ

    @staticmethod
    def _parse_logfile(logcontents):
        """ Method to parse the logcontents and return eLLG and VRMS

        :param logcontents: the contents of the log file
        :type logcontents: str
        :returns a tuple with the eLLG (str) and the VRMS (str)
        :rtype: tuple
        """

        eLLG = "NA"
        VRMS = "NA"
        # Parse the logfile for eLLG and VRMS
        ellg_reached = False
        for line in logcontents.split("\n"):
            line = line.rstrip().lstrip()
            if "eLLG   RMSD frac-scat  Ensemble" in line:
                ellg_reached = True
            elif ellg_reached:
                eLLG = line.split()[0]
                ellg_reached = False
            if "SOLU ENSEMBLE" in line and "VRMS DELTA" in line:
                VRMS = line.split()[5].rstrip().lstrip()
                break
        return eLLG, VRMS
