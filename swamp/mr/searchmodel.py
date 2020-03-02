import abc
import os
import gemmi
import shutil
import logging
import swamp
from swamp.utils import decompress
from swamp.wrappers import Gesamt
from swamp.utils import renumber_hierarchy

ABC = abc.ABCMeta('ABC', (object,), {})


class SearchModel(ABC):
    """Class with methods to prepare the search model before MR and data structures to keep all useful information

    :param str workdir: working directory where the search model will be prepared
    :param str id: unique identifier for the search model to be added
    :param str ensemble_code: the ensemble's SWAMP library id to be used as search model
    :param float ermsd: the eRMSD to be used with phaser to place the search model (default 0.1)
    :param int nsearch: number of copies to search with phaser
    :param bool disable_check: passed to :py:obj:`phaser.InputMR_AUTO.setENSE_DISA_CHEC` (default True)
    :param str mod: indicate how to prepare the search model (default 'unmod')
    :param str model: indicate if the search model is an ensemble or a centroid (default 'ensemble')
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    :ivar bool error: if True an error has occurred while preparing the search model
    :ivar list modified_model_list: a list with the file names of the modified models to be merged into the new ensemble
    """

    def __init__(self, id, ensemble_code, workdir, ermsd=0.1, nsearch=1, disable_check=True, mod='unmod',
                 model='ensemble', logger=None):

        self.id = id
        self.ensemble_code = ensemble_code
        self.workdir = workdir
        self.ermsd = ermsd
        self.disable_check = disable_check
        self.nsearch = nsearch
        self.mod = mod
        self.model = model
        self.error = False
        self.modified_model_list = []

        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        if ensemble_code == 'idealhelix':
            self.gzfile = os.path.join(self.idealhelix_fname)
            self._pdbfname = os.path.join(self.workdir, 'idealhelix.pdb')
        else:
            self.gzfile = os.path.join(swamp.ENSEMBLE_DIR, '%s_%s.pdb.gz' % (model, ensemble_code))
            self._pdbfname = os.path.join(self.workdir, '%s_%s.pdb' % (model, ensemble_code))

        if not os.path.isfile(self.gzfile):
            self.logger.error('Search model file not found! %s\nMake sure the ensemble code is correct!' % self.gzfile)
            self.error = True
            return

        self._make_workdir()
        self.create_pdbfile()

    # ------------------ Properties ------------------

    @property
    def phaser_info(self):
        """A dictionary with all the information necessary to use with \
        :py:func:`~swamp.wrappers.wphaser.Phaser.add_searchmodel`"""

        return {'id': self.id,
                'pdbfile': self.pdbfname,
                'ermsd': self.ermsd,
                'nsearch': self.nsearch,
                'disable_check': self.disable_check}

    @property
    def idealhelix_fname(self):
        """File name of the ideal helix to be used to extend the solution"""
        return os.path.join(swamp.IDEALHELICES_DIR, 'ensemble_20_nativebfact_homogenous.pdb.gz')

    @property
    def pdbfname(self):
        """The PDB file name of the search model"""
        if self.mod == 'unmod' or self.ensemble_code == 'idealhelix':
            return self._pdbfname
        else:
            return self.modified_pdbfname

    @property
    def model_list(self):
        """A list with all the models in the original search model"""
        if not os.path.isdir(self.model_dir):
            os.makedirs(self.model_dir)
        return self.split_models(pdbin=self._pdbfname, directory=self.model_dir, strip_hetatm=True)

    @property
    def model_dir(self):
        """Directory with the models that formed the original search model (if it was an ensemble)"""
        return os.path.join(self.workdir, "models")

    @property
    def modified_pdbfname(self):
        """The file name of the search model after running :py:func:`~swamp.mr.searchmodel.prepare`"""
        return os.path.join(self.workdir, '%s_%s_%s.pdb' % (self.mod, self.model, self.ensemble_code))

    @property
    def _modified_model_template(self):
        """String to be used as a template for the modified model file name"""
        return os.path.join(self.workdir, "{}_%s.pdb" % self.mod)

    # ------------------ Methods ------------------

    def create_pdbfile(self):
        """Create the pdb file to be used as a search model on :py:attr:`~swamp.mr.mrrun.MrRun.run`"""

        decompress(self.gzfile, self._pdbfname)
        if self.mod == 'unmod' or self.ensemble_code == 'idealhelix':
            return
        self.prepare()
        self._check_output()

    def prepare(self):
        """Prepare the :py:attr:`~swmap.mr.searchmodel.Searchmodel.pdbfname` with the searchmodel using the indicated \
        :py:attr:`~swmap.mr.searchmodel.Searchmodel.mod` before MR"""

        if self.mod == 'polyala':
            for idx, model in enumerate(self.model_list):
                modelID = os.path.basename(model)[:-4]
                modified_model = self._modified_model_template.format(modelID)
                self.logger.debug('Truncating model %s %s -> %s' % (idx, modelID, modified_model))
                self.truncate_polyALA(pdbin=model, pdbout=modified_model)
                self.logger.debug('Transfer flags to new pdb file')
                self.transfer_flags_pdb(pdb_ref=model, pdb_file=modified_model)
                self.modified_model_list.append(modified_model)
            self._merge_models()

        else:
            self.extract_core(pdbout=self.modified_pdbfname, workdir=self.workdir, model_list=self.model_list)

    # ------------------ Hidden methods ------------------

    def _check_output(self):
        """Check if the output file has been created, set :py:attr:`~swamp.searchmodel.prepare.error` to \
        True if not"""

        if self.mod != 'unmod' and not os.path.isfile(self.modified_pdbfname):
            self.logger.error("Modified search model not found! %s" % self.modified_pdbfname)
            self.error = True

    def _merge_models(self):
        """Merge all the modified models indicated at \
        :py:attr:`~swamp.searchmodel.prepare.modified_model_list` into \
        :py:attr:`~swamp.searchmodel.prepare.pdbout`"""

        if len(self.model_list) > 1:
            gesamt = Gesamt(mode='alignment', pdbin=self.modified_model_list, pdbout=self.modified_pdbfname,
                            workdir=self.workdir, logger=self.logger)
            gesamt.run()

        else:
            shutil.copyfile(self.modified_model_list[0], self.modified_pdbfname)

    def _make_workdir(self):
        """Method to crete the :py:attr:`~swamp.searchmodel.prepare.workdir` directory"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    # ------------------ Static methods ------------------

    @staticmethod
    def extract_core(pdbout, workdir, model_list):
        """Extract the core structural alignment of a given ensemble

        :param str pdbout: the pdb file name of the resulting ensemble to output
        :param str workdir: the working directory where temporary files will be created
        :param list model_list: a list with the pdb file names of the models that make the ensemble
        """

        from mrbump.seq_align.MRBUMP_gesamt import Gesamt as MRBUMP_Gesamt

        tmp_csvfile = os.path.join(workdir, "csvfile.csv")
        tmp_csvCOREfile = os.path.join(workdir, "csvCOREfile.csv")
        tmp_logfile = os.path.join(workdir, "core_logfile.log")
        tmp_alnfile = os.path.join(workdir, "alnfile.ali")
        tmp_scriptfile = os.path.join(workdir, "scriptfile.sh")
        tmp_pdb = os.path.join(workdir, "tmp_out_core.pdb")

        pdbDict = dict([])
        for model in model_list:
            if len(os.path.basename(model)[:-4]) >= 6:
                pdbDict[model] = "%s" % os.path.basename(model)[:-4][-1]
            else:
                pdbDict[model] = "*"

        ensemble_truncator = MRBUMP_Gesamt()
        ensemble_truncator.runGesamt(model_list, pdbDict, seqin="fake.fasta", outputPDB=tmp_pdb,
                                     logfile=tmp_logfile, alnfile=tmp_alnfile, csvfile=tmp_csvfile,
                                     script=tmp_scriptfile, debug=False)

        ensemble_truncator.makeGesTruncEnsemble(tmp_pdb, pdbout, variancePercent=10000, sidechain_level=100,
                                                csvFile=tmp_csvfile, truncation_level=1000.0,
                                                csvCOREfile=tmp_csvCOREfile)

    @staticmethod
    def split_models(pdbin, directory, strip_hetatm=True):
        """Method to split an ensemble into its model components

        :param str pdbin: input pdb file name
        :param str directory: directory where the models of the ensemble will be dumped
        :param bool strip_hetatm: if set, the hetatm will be ommited from the output models (default True)
        :returns: a list with the output file names listed (list)
        """

        pdbout_template = os.path.join(directory, 'model_{}.pdb')
        hierarchy = gemmi.read_structure(pdbin)
        rslt = []

        if strip_hetatm:
            hierarchy.remove_waters()

        for model in hierarchy:
            new_structure = gemmi.Structure()
            new_structure.add_model(model)
            new_structure.cell = hierarchy.cell
            new_structure.write_pdb(pdbout_template.format(model.name))
            rslt.append(pdbout_template.format(model.name))

        return rslt

    @staticmethod
    def truncate_polyALA(pdbin, pdbout):
        """Method to truncate a given pdb into poly-ala

        :param str pdbin: input pdb file name
        :param str pdbout: output pdb file name
        """

        original_hierarchy = gemmi.read_structure(pdbin)
        original_hierarchy.remove_ligands_and_waters()

        for residue in original_hierarchy[0][0]:
            residue.trim_to_alanine()
            residue.name = "ALA"

        renumber_hierarchy(original_hierarchy)

        original_hierarchy.write_pdb(pdbout)

    @staticmethod
    def transfer_flags_pdb(pdb_ref, pdb_file, flags_to_transfer=("CRYST1", "SCALE", "REMARK"), overwrite=True):
        """Transfer PDB flags between two given pdb files

        :param str pdb_ref: pdb file with the reference flags to be transferred
        :param str pdb_file: pdb file where the flags will be transferred
        :param tuple flags_to_transfer: set of flags that need to be transferred
        :param bool overwrite: if False, pdb_file original flags will be  kept (default True)
        """

        with open(pdb_ref, "r") as pdbreference, open(pdb_file, "r") as pdbfile:
            # Read in the lines to transfer
            lines_to_include = []
            for line in pdbreference:
                if line.split()[0] in flags_to_transfer:
                    lines_to_include.append(line)
            # Store the lines of the pdbfile
            for line in pdbfile:
                if not overwrite:
                    lines_to_include.append(line)
                elif line.split()[0] not in flags_to_transfer:
                    lines_to_include.append(line)
        # Re-open the file, this time write mode
        pdbreference.close()
        pdbfile.close()
        with open(pdb_file, "w") as pdbfile:
            for line in lines_to_include:
                pdbfile.write(line)
