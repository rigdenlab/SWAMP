import os
import shutil
import gemmi
import shlex
import subprocess
from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.SeqUtils import molecular_weight
from simbad.util.mtz_util import GetLabels
from simbad.parsers.anode_parser import AnodeParser
from simbad.mr.anomalous_util import AnodeSearch
from mmtbx.scaling.matthews import matthews_rupp
from cctbx.crystal import symmetry
from iotbx import reflection_file_reader
import swamp
from swamp.mr.mr import Mr
from swamp.wrappers.shelxe import Shelxe
from swamp.wrappers.crank2 import Crank2
from swamp.wrappers.wphaser import Phaser
from swamp.wrappers.wrefmac import wRefmac
from swamp.library.tools import decompress
from swamp.searchmodel_prepare.core import Core
from swamp.searchmodel_prepare.polyala import PolyALA
from swamp.searchmodel_prepare.bfactor import Bfactor
from swamp.searchmodel_prepare.molrep import Molrep
from swamp.searchmodel_prepare.centroid import Centroid


class MrRun(Mr):
    """Class to attempt molecular replacement on a given target

    This class contains methods to attempt run molecular replacement on a given target structure. The pipeline uses
    phaser for search model placement, refmac5 for refinement, shelxe for density modification and model building. If
    scattering data has been recorded into the target's mtz file, anode is used to detect anomalous peaks followed by
    crank2 for experimental phasing.

    :param str, int id: unique identifier for the instance
    :param str workdir: working directory where the MR task will be executed
    :param int threads: number of threads to be used in the pipeline (only affects phaser) (default 1)
    :param str target_fa: target's fasta filename
    :param str target_mtz: target's mtz filename
    :param str phased_mtz: target's mtz filename containing phases (default: None)
    :param str phaser_sgalternative: parameter to be passed to phaser as SGAL_SELE (default 'NONE')
    :param float, None phaser_packcutoff: parameter to be passed to phaser as setPACK_CUTO (default None)
    :param float, None phaser_peaks_rotcutoff: parameter to be passed to phaser as setPEAK_ROTA_CUTO (default None)
    :param bool phaser_early_kill: pipeline will stop execution if phaser placement is incorrect when True (default False)
    :param int phaser_timeout: parameter to be passed to phaser as KILL_TIME (default 360)
    :param bool extend_solution: if True then solutions will be completed with ideal helices (default False)
    :param bool save_disk_space: if True metadata will be removed to save disk space (default False)
    :param bool ignore_anomalous: if True then presence of scattering data will not trigger anode/crank2 execution
    :param logger: logging interface for the MR pipeline
    :param bool silent: if set to True the logger will not print messages
    :param bool quiet_start: if True the logger will not display the header section and the inital parameters
    :ivar error: True if errors have occurred at some point on the pipeline
    :ivar is_extended: indicates number of cycles used to extend the solution
    :ivar phaser: phaser wrapper used in the pipeline
    :ivar refmac: refmac wrapper used in the pipeline
    :ivar shelxe: shelxe wrapper used in the pipeline
    :ivar anode: anode wrapper used in the pipeline
    :ivar crank2: crank2 wrapper used in the pipeline
    :ivar anode_parser: anode lsa output file parser wrapper
    """

    def __init__(self, id, workdir, target_fa, target_mtz, phased_mtz=None, threads=1, phaser_sgalternative="NONE",
                 phaser_early_kill=True, ignore_anomalous=True, silent=False, save_disk_space=False, logger=None,
                 phaser_timeout=1440, extend_solution=True, quiet_start=False, phaser_packcutoff=None,
                 phaser_peaks_rotcutoff=None):

        super(MrRun, self).__init__(id, target_fa, target_mtz, workdir, phased_mtz=phased_mtz, logger=logger,
                                    silent=silent)

        self.init_params = locals()
        if not quiet_start:
            self.logger.info(self.pipeline_header.format(' MR-RUN '))
            self.logger.info(self._inform_args(**self.init_params))
        self._searchmodel_list = []
        self._threads = threads
        self._save_disk_space = save_disk_space
        self._extend_solution = extend_solution
        self._phaser_sgalternative = phaser_sgalternative
        self._phaser_early_kill = phaser_early_kill
        self._phaser_timeout = phaser_timeout
        self._ignore_anomalous = ignore_anomalous
        self._is_extended = "NO"
        self._phaser = None
        self._refmac = None
        self._shelxe = None
        self._crank2 = None
        self._anode = None
        self._anode_parser = None
        self._parent_array = None
        self._solution = None
        self._idealhelix_run = None
        self._phaser_packcutoff = phaser_packcutoff
        self._phaser_peaks_rotcutoff = phaser_peaks_rotcutoff

    # ------------------ Some general properties ------------------

    @property
    def cleanup_dir_list(self):
        """Property storing information about directories used in the pipeline that can be removed to save disk space"""

        return [self.phaser_info['workdir'], self.refmac_info['workdir'], self.shelxe_info['workdir'],
                os.path.join(self.workdir, "anode"), self.crank2_info['workdir'], self.searchmodel_dir,
                os.path.join(self.workdir, "ideal_helices")]

    @property
    def save_disk_space(self):
        """True if it is necessary to save disk space and remove unnecessary results"""
        return self._save_disk_space

    @save_disk_space.setter
    def save_disk_space(self, value):
        self._save_disk_space = value

    @property
    def idealhelix_run(self):
        """Instance of :obj:`swamp.mr.mrrun.MrRun to extend the solution with ideal helices"""
        return self._idealhelix_run

    @idealhelix_run.setter
    def idealhelix_run(self, value):
        self._idealhelix_run = value

    @property
    def phaser_packcutoff(self):
        return self._phaser_packcutoff

    @phaser_packcutoff.setter
    def phaser_packcutoff(self, value):
        self._phaser_packcutoff = value

    @property
    def phaser_peaks_rotcutoff(self):
        return self._phaser_peaks_rotcutoff

    @phaser_peaks_rotcutoff.setter
    def phaser_peaks_rotcutoff(self, value):
        self._phaser_peaks_rotcutoff = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        self._threads = value

    @property
    def phaser_early_kill(self):
        return self._phaser_early_kill

    @phaser_early_kill.setter
    def phaser_early_kill(self, value):
        self._phaser_early_kill = value

    @property
    def phaser_sgalternative(self):
        return self._phaser_sgalternative

    @phaser_sgalternative.setter
    def phaser_sgalternative(self, value):
        self._phaser_sgalternative = value

    @property
    def phaser_timeout(self):
        return self._phaser_timeout

    @phaser_timeout.setter
    def phaser_timeout(self, value):
        self._phaser_timeout = value

    @property
    def phaser(self):
        """Property to hold the phaser wrapper :obj:`swamp.wrappers.wphaser.Phaser` instance"""
        return self._phaser

    @phaser.setter
    def phaser(self, value):
        self._phaser = value

    @property
    def crank2(self):
        """Property to hold the crank2 wrapper :obj:`swamp.wrappers.crank2.Crank2` instance"""
        return self._crank2

    @crank2.setter
    def crank2(self, value):
        self._crank2 = value

    @property
    def anode(self):
        """Property to hold the anode wrapper"""
        return self._anode

    @anode.setter
    def anode(self, value):
        self._anode = value

    @property
    def shelxe(self):
        """Property to hold the shelxe wrapper :obj:`swamp.wrappers.shelxe.Shelxe` instance"""
        return self._shelxe

    @shelxe.setter
    def shelxe(self, value):
        self._shelxe = value

    @property
    def refmac(self):
        """Property to hold the refmac wrapper :obj:`swamp.wrappers.wrefmac.Refmac` instance"""
        return self._refmac

    @refmac.setter
    def refmac(self, value):
        self._refmac = value

    @property
    def extend_solution(self):
        """If True extend the solution with ideal helices"""
        return self._extend_solution

    @extend_solution.setter
    def extend_solution(self, value):
        self._extend_solution = value

    @property
    def ignore_anomalous(self):
        """If True presence of anomalous signal must be ignored"""
        return self._ignore_anomalous

    @ignore_anomalous.setter
    def ignore_anomalous(self, value):
        self._ignore_anomalous = value

    @property
    def searchmodel_list(self):
        """List of the search models to be used in this Mr run"""
        return self._searchmodel_list

    @searchmodel_list.setter
    def searchmodel_list(self, value):
        self._searchmodel_list = value

    @property
    def is_extended(self):
        """If True ideal helices have been used to extend the solution"""
        return self._is_extended

    @is_extended.setter
    def is_extended(self, value):
        self._is_extended = value

    @property
    def solution(self):
        """If 'YES' the MrRun solved the structure"""
        return self._solution

    @solution.setter
    def solution(self, value):
        self._solution = value

    @property
    def anode_parser(self):
        """Property to hold the anode parser"""
        return self._anode_parser

    @anode_parser.setter
    def anode_parser(self, value):
        self._anode_parser = value

    @property
    def searchmodel_dir(self):
        """Directory where the search model preparation will take place"""
        return os.path.join(self.workdir, "searchmodels")

    @property
    def prepared_searchmodel_fname(self):
        """Filename of the modified search model"""
        return os.path.join(self.searchmodel_dir, "searchmodel_{}_{}.pdb")

    @property
    def target_info(self):
        """Dictionary to contain information about the target

        :raises ValueError if the target's mtz and fasta files were not provided
        :returns a dictionary with the target information
        :rtype dict
        """

        if self.target_mtz is None or self.target_fa is None:
            self.error = True
            raise ValueError('Need to setup target mtz and fasta files before running MR pipeline!')

        MW, solvent, nchains_asu, anomalous_signal, nresidues_chain, \
        resolution, spacegroup, spacegroup_symbol, n_reflections, \
        use_f = self.characterise_target(self.target_fa, self.target_mtz)

        return {'fastafile': self.target_fa,
                'mtzfile': self.target_mtz,
                'resolution': resolution,
                'use_f': use_f,
                'MW': MW,
                'phasedmtz': self.phased_mtz,
                'solvent': solvent,
                'nchains_asu': nchains_asu,
                'anomalous_signal': anomalous_signal,
                'nresidues_chain': nresidues_chain,
                'n_reflections': n_reflections,
                }

    @property
    def phaser_info(self):
        """Dictionary to contain information about the phaser run

        :returns a dictionary with the arguments for phaser
        :rtype dict
        """

        return {'early_kill': self.phaser_early_kill,
                'workdir': os.path.join(self.workdir, "phaser"),
                'timeout': self.phaser_timeout,
                'logger': self.logger,
                'threads': self.threads,
                'phased_mtz': self.phased_mtz,
                'mtzfile': self.target_mtz,
                'mw': self.target_info['MW'],
                'packcutoff': self.phaser_packcutoff,
                'nchains_asu': self.target_info['nchains_asu'],
                'sgalternative': self.phaser_sgalternative,
                'peaks_rotcutoff': self.phaser_peaks_rotcutoff
                }

    @property
    def refmac_info(self):
        """Dictionary to contain information about the refmac run

        :returns a dictionary with the arguments for refmac
        :rtype dict
        """

        return {'workdir': os.path.join(self.workdir, "refmac"),
                'pdbin': self.phaser.pdbout,
                'mtzin': self.target_mtz,
                'phased_mtz': self.phased_mtz,
                'logger': self.logger
                }

    @property
    def crank2_info(self):
        """Dictionary to contain information about the crank2 run

        :returns a dictionary with the arguments for crank2
        :rtype dict
        """

        return {'workdir': os.path.join(self.workdir, "crank2"),
                'logger': self.logger,
                'pdbin': os.path.join(self.refmac.pdbout),
                'mtzin': self.target_mtz,
                'anomalous_scatterer': "S",
                'mode': "SAD",
                'wavelength_type': "peak",
                'fastafile': self.target_fa,
                'nchains_asu': self.target_info['nchains_asu'],
                'solvent': self.target_info['solvent']
                }

    @property
    def shelxe_info(self):
        """Dictionary to contain information about the shelxe run

        :returns a dictionary with the arguments for shelxe
        :rtype dict
        """

        return {'workdir': os.path.join(self.workdir, 'shelxe'),
                'logger': self.logger,
                'pdbin': self.refmac.pdbout,
                'mtzin': self.target_mtz,
                'solvent': self.target_info['solvent'],
                'nreflections': self.target_info['n_reflections'],
                'use_f': self.target_info['use_f'],
                'resolution': self.target_info['resolution']
                }

    @property
    def _list_idealhelices(self):
        """List of the ideal helices available to extend the solution"""

        permited_sizes = ["10", "15", "20", "25"]
        permited_modification = ["nativebfact", "gradientbfact"]
        return [x for x in os.listdir(swamp.IDEALHELICES_DIR) if
                x.split("_")[1] in permited_sizes
                and x.split("_")[-2] in permited_modification
                and x.split("_")[-1] == "homogenous.pdb"]

    @property
    def idealhelix_fname(self):
        """File name of the ideal helix to be used to extend the solution"""
        return os.path.join(swamp.IDEALHELICES_DIR, 'ensemble_20_nativebfact_homogenous.pdb.gz')

    @property
    def idealhelices_workdir(self):
        """Directory where the ideal helices solution extension will take place"""
        return os.path.join(self.workdir, "ideal_helices")

    # ------------------ Some general methods ------------------

    def add_searchmodel(self, id, ensemble_code, ermsd=0.1, nsearch=1, disable_check=True, modification='unmod'):
        """Add a search model to the phaser run

        :param id: unique identifier for the search model to be added
        :type id: str, int
        :param ensemble_code: the ensemble's SWAMP library id to be used as search model
        :type ensemble_code: str, int
        :param ermsd: the eRMSD to be used with phaser to place the search model (default 0.1)
        :type ermsd: float
        :param nsearch: number of copies to search with phaser
        :type nsearch: int
        :param disable_check: passed to :obj:`phaser.InputMR_AUTO.setENSE_DISA_CHEC` (default True)
        :type disable_check: bool
        :param modification: indicate how to prepare the search model (default 'unmod')
        :type modification: str
        :returns nothing
        :rtype None
        """

        if ensemble_code == 'idealhelix':
            gzfile = os.path.join(self.idealhelix_fname)
            pdbfile = os.path.join(self.workdir, 'idealhelix.pdb')
        else:
            gzfile = os.path.join(swamp.ENSEMBLE_DIR, 'ensemble_%s.pdb.gz' % ensemble_code)
            pdbfile = os.path.join(self.workdir, 'ensemble_%s.pdb' % ensemble_code)

        if not os.path.isfile(gzfile):
            self.logger.error('Search model file not found! %s\nMake sure the ensemble code is correct!' % gzfile)
            self.error = True
            return

        decompress(gzfile, pdbfile)

        if self.searchmodel_list and id in [x['id'] for x in self.searchmodel_list]:
            self.logger.error('A searchmodel with the same id has been already added!')
            self.error = True
            return

        if modification != 'unmod':
            fname = self.prepared_searchmodel_fname.format(id, modification)
            self.prepare_search_model(modification=modification, pdbin=pdbfile, workdir=self.searchmodel_dir,
                                      pdbout=fname)
        else:
            fname = pdbfile

        if not self.error:
            self.searchmodel_list.append({
                'id': id,
                'pdbfile': fname,
                'ermsd': ermsd,
                'nsearch': nsearch,
                'disable_check': disable_check
            })
        else:
            self.logger.warning('Previous errors prevent adding the searchmodel!')
            return

    def register_solution(self, **kwargs):
        """Register an existing solution information to be used with phaser"""
        self.solution = kwargs

    def append_results(self):
        """Method to append the results obtained into the result list"""

        self.results.append(
            [self.id, self.phaser.LLG, self.phaser.TFZ, self.phaser.local_CC, self.phaser.overall_CC,
             self.refmac.rfree, self.refmac.rfactor, self.refmac.local_CC, self.refmac.overall_CC, self.shelxe.cc,
             self.shelxe.acl, self.anode.peak_height, self.crank2.rfree, self.crank2.rfactor, self.is_extended,
             self.shelxe.solution])

    def run(self):
        """Run the MR pipeline using of phaser, refmac5 and shelxe. Extend the possible solution with ideal helices.

        :returns no value
        :rtype None
        """

        self._initiate_wrappers()

        # Sanity check
        if self.error:
            self.logger.warning("Previous errors prevent execution of the pipeline")
            return

        # Run phaser
        self.phaser.run()

        # If there is a problem, abort
        if self.phaser.error or self.phaser.abort_suggested:
            if self.phaser.error:
                self.logger.warning("Previous error prevents pipeline moving forward...")
            else:
                self.logger.warning("Phaser scores below threshold, early termination triggered...")
            self.append_results()
            return

        # Run refmac
        self.refmac.run()
        self.refmac.make_logfile()

        # If didn't give output, abort
        if self.refmac.error:
            self.logger.warning("Previous error prevents pipeline moving forward...")
            self.append_results()
            return

        # If there is anomalous signal, look for it
        if self.target_info['anomalous_signal'] and not self.ignore_anomalous and self._found_anomalous_signal():
            self.crank2.run()
            if self.crank2.solution == "YES":
                self.logger.info("Crank2 found a solution!")
                self.append_results()
                return

        # Run shelxe
        self.shelxe.run()
        self.shelxe.make_logfile()

        # Store the results
        self.append_results()

        # If there is no solution and the cc is promising, let's try to fit some individual helices
        if self.extend_solution and self.shelxe.solution == "NO" and self.crank2.solution == "NO":
            self.logger.info("Search model placement did not yield a solution, trying to fit ideal helices!\n")
            self.fit_helices()

        # If there was no solution and the user wants to save disk space, delete unnecessary stuff
        if self.shelxe.solution == "NO" and self.save_disk_space and self.idealhelix_run.shelxe.solution == 'NO':
            self.logger.info("Saving disk space, %s will be deleted!" % self.workdir)
            os.chdir(os.path.dirname(self.workdir))
            self._cleanup_files()

        # Exit and show the results unless this was an ideal helix run
        if not self.is_extended:
            self.logger.info('MR run finished. Table of results:\n\n%s\n' % self.table_contents)

    def fit_helices(self):
        """Method to extend the solution with ideal helices

        This method will create and set running another :obj:`swamp.mr.mrrun.MrRun` instance that will take
        the placed search model as an existing solution and try to place ideal helices to extend it.
        """

        # Set up input parameters based on original MR run parameters with some tweaks
        size = 20
        n_helices = int(self.target_info['nresidues_chain'] / size)
        n_helices -= 2 * len(self.searchmodel_list)
        input_params = {'id': self.id,
                        'workdir': os.path.join(self.idealhelices_workdir),
                        'target_fa': self.target_fa,
                        'target_mtz': self.target_mtz,
                        'phased_mtz': self.phased_mtz,
                        'phaser_sgalternative': self.phaser_sgalternative,
                        'phaser_packcutoff': 40.0,
                        'phaser_timeout': 1440,
                        'silent': False,
                        'ignore_anomalous': True,
                        'phaser_early_kill': True,
                        'logger': self.logger,
                        'extend_solution': False,
                        'save_disk_space': False,
                        'quiet_start': True,
                        'threads': self.threads
                        }
        # Run  MR pipeline using the ideal helix as search model and manage results
        self.idealhelix_run = MrRun(**input_params)
        self.idealhelix_run.add_searchmodel('idealhelix', ensemble_code='idealhelix', nsearch=5)
        self.idealhelix_run.is_extended = 'YES'
        self.idealhelix_run.register_solution(pdbfile=self.refmac.pdbout, ermsd=0.1)
        self.idealhelix_run.run()
        self.results += self.idealhelix_run.results

    def prepare_search_model(self, modification='polyALA', **kwargs):
        """ Method to prepare the search model with a given modification protocol

        :param modification: indicates the modification to be used (default 'polyALA')
        :type modification: str
        :param kwargs: arguments to be passed to :obj:`swamp.searchmodel_prepare.prepare`
        :type kwargs: dict
        :returns nothing
        :rtype None
        """

        self.logger.info('Preparing searchmodel pdbfile for MR run')

        if modification == "unmod":
            shutil.copyfile(kwargs['pdbin'], kwargs['pdbout'])
            return
        elif modification == "polyALA":
            prep = PolyALA(**kwargs)
        elif modification == "core":
            prep = Core(**kwargs)
        elif modification == "bfactor":
            prep = Bfactor(**kwargs)
        elif modification == "centroid":
            prep = Centroid(**kwargs)
        elif modification == "molrep":
            prep = Molrep(**kwargs)

        prep.prepare()

        if prep.error:
            self.error = True
            self.logger.warning("Search model modification wasn't completed due to previous errors...")

            print(modification)
            print(self.id)
            print(self.searchmodel_list)

    def characterise_searchmodel(self, searchmodel_filename):
        """Method to characterise a given search model (no. of models, no. of residues and avg. qscore of the ensemble)
    
        :param searchmodel_filename: filename of the pdb file of the search model
        :type searchmodel_filename: str
        :returns tuple containing the number of models, number of residues and average qscore of the ensemble
        :rtype: tuple
        """

        # Read hierarchy
        hierarchy = gemmi.read_structure(searchmodel_filename)
        hierarchy.remove_ligands_and_waters()

        # Get the qscore
        gesamtEXE = os.path.join(os.environ["CCP4"], "bin", "gesamt")
        for model in hierarchy:
            gesamtEXE += " %s -s /%s/" % (searchmodel_filename, model.name)
        process_args = shlex.split(gesamtEXE, posix=False)
        p = subprocess.Popen(process_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
        avg_qscore = "NA"
        for line in stdout.split("\n"):
            if "Q-score" in line:
                avg_qscore = line.rstrip().lstrip().split(":")[1].split()[0].rstrip().lstrip()
                break
        if avg_qscore == "NA":
            self.logger.warning("Something went wrong getting the avg. Qscore of the search model!")

        # Return data
        return len(hierarchy), len([residue for model in hierarchy for chain in model for residue in chain]), avg_qscore

    # ------------------ Some hidden methods ------------------

    def _initiate_wrappers(self):
        """Method to instantiate the wrappers to be used in the pipeline

        This method will instantiate the :onj:`swamp.wrapper.wrapper` classes required for the pipeline execution:
        phaser, refmac, shelxe, crank2, anode, phenix_get_cc
        """

        if not self.searchmodel_list:
            self.error = True
            self.logger.error('Cannot proceed with MR run without at least one search model!')
            return

        self.phaser = Phaser(**self.phaser_info)
        for searchmodel_info in self.searchmodel_list:
            self.phaser.add_searchmodel(**searchmodel_info)
        if self.solution is not None:
            self.phaser.register_solution(**self.solution)
        self.refmac = wRefmac(**self.refmac_info)
        self.shelxe = Shelxe(**self.shelxe_info)
        self.crank2 = Crank2(**self.crank2_info)
        if not os.path.isdir(os.path.join(self.workdir, "anode")) and not self.ignore_anomalous:
            os.mkdir(os.path.join(self.workdir, "anode"))
        self._anode = AnodeSearch(os.path.abspath(self.target_info['mtzfile']), os.path.join(self.workdir, "anode"))
        self.anode.peak_height = "NA"

    def _found_anomalous_signal(self):
        """Look for the presence of an anomlaous peak using anode"""

        self.logger.info("Looking for anomalous signal using ANODE")
        self.anode.run(self.refmac.pdbout)
        self.anode.peak_height = self.get_anomalous_peak(
            os.path.join(self.workdir, "anode", "%s.lsa" % os.path.basename(self.refmac.pdbout)[:-4]))
        self.logger.info("Anomalous peak height: %s" % self.anode_parser.peak_height)
        if self.anode.peak_height != "NA" and float(self.anode.peak_height) > 7.29:
            return True
        else:
            return False

    # ------------------ Some static methods ------------------

    @staticmethod
    def get_anomalous_peak(lsa_file):
        """Method to detect anomalous signal peak in a lsa file
    
        :param lsa_file: filename of the lsa anode output to be parsed
        :type lsa_file: str
        :returns peak_height: intensity of the highest peak detected ('NA' if there are no peaks)
        :rtype str
        """

        anode_parser = AnodeParser(lsa_file)
        if anode_parser.peak_height is None:
            return 'NA'
        else:
            return str(anode_parser.peak_height)

    @staticmethod
    def characterise_target(target_fa, target_mtz):
        """Method to characterise a given target
    
        :param target_fa: filename of the fasta file of the target
        :type target_fa: str
        :param target_mtz: filename of the mtz file of the target
        :type target_mtz: str
        :returns tuple containing target's characteristics
        :rtype: tuple
        """

        # Sequence information
        target_chains = [str(chain.seq) for chain in list(SeqIO.parse(target_fa, "fasta", alphabet=generic_protein))]
        target_chains = list(set(target_chains))
        target_molecular_weight = 0.0
        for seq in target_chains:
            seq = seq.replace("X", "A")
            target_molecular_weight += round(molecular_weight(seq, "protein"), 2)
        seq_length = sum([len(seq) for seq in target_chains])

        # Crystal information
        mtz_head = GetLabels(target_mtz)
        if mtz_head.i is None and mtz_head.f is not None:
            use_f = True
        else:
            use_f = False
        reflection_file = reflection_file_reader.any_reflection_file(file_name=target_mtz)
        content = reflection_file.file_content()
        n_reflections = int(content.n_reflections())
        target_space_group = content.space_group_name().replace(' ', '')
        spacegroup_symbol = content.space_group_info().symbol_and_number()
        spacegroup_symbol = spacegroup_symbol.split("(No.")[0].rstrip().lstrip()
        if ":" in spacegroup_symbol:
            spacegroup_symbol = spacegroup_symbol.split(":")[0].rstrip().lstrip()
        resolution = float(content.max_min_resolution()[1])
        target_unit_cell = content.crystals()[0].unit_cell_parameters()
        crystal_symmetry = symmetry(space_group_symbol=target_space_group, unit_cell=target_unit_cell)
        matthews_coeff = matthews_rupp(crystal_symmetry, n_residues=seq_length)
        solvent = round(matthews_coeff.solvent_content, 2)
        ncopies = matthews_coeff.n_copies
        mtz_head = GetLabels(target_mtz)
        if mtz_head.sigiplus is not None and mtz_head.iplus is not None and mtz_head.sigiminus is not None \
                and mtz_head.iminus is not None:
            anomalous_signal = True
        else:
            anomalous_signal = False

        return target_molecular_weight, solvent, ncopies, anomalous_signal, seq_length, resolution, \
               target_space_group, spacegroup_symbol, n_reflections, use_f
