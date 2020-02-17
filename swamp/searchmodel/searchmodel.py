import abc
import os
import gemmi
import shutil
import logging
from swamp.wrappers import Gesamt

ABC = abc.ABCMeta('ABC', (object,), {})


class SearchModel(ABC):
    """Abstract class for serch model preparation

    :param str workdir: working directory where the search model will be prepared
    :param str pdbin: input pdb file name
    :param str pdbout: output pdb file name
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    :ivar bool error: if True an error has occurred while preparing the search model
    :ivar str logcontents: stores the stdout of any command line tool executed while preparing the search model
    :ivar list modified_model_list: a list with the file names of the modified models to be merged into the pdbout
    """

    def __init__(self, workdir, pdbin, pdbout, logger=None):
        self.workdir = workdir
        self.pdbin = pdbin
        self.pdbout = pdbout
        self.error = False
        self.logcontents = None
        self.modified_model_list = []
        if logger is None:
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def prepare(self):
        """Abstract method to prepare the search model"""
        pass

    @property
    @abc.abstractmethod
    def modification(self):
        """Abstract property to store modification to be applied"""
        pass

    # ------------------ Prroperties ------------------

    @property
    def model_list(self):
        """A list with all the models in the original search model"""
        if not os.path.isdir(self.model_dir):
            os.makedirs(self.model_dir)
            return self.split_models(pdbin=self.pdbin, directory=self.model_dir, strip_hetatm=True)
        else:
            return [os.path.join(self.model_dir, filename) for filename in os.listdir(self.model_dir)]

    @property
    def model_dir(self):
        """Directory with the models that formed the original search model (if it was an ensemble)"""
        return os.path.join(self.workdir, "models")

    @property
    def _modified_model_template(self):
        """String to be used as a template for the modified model file name"""
        return os.path.join(self.workdir, "{}_%s.pdb" % self.modification)

    @property
    def _tmp_pdb(self):
        """Temporary pdb file name"""
        return os.path.join(self.workdir, "tmp_out_%s.pdb" % self.modification)

    # ------------------ Hidden methods ------------------

    def _check_output(self):
        """Check if the output file has been created, set :py:attr:`~swamp.searchmodel.prepare.error` to \
        True if not"""

        if not os.path.isfile(self.pdbout):
            self.logger.error("Modified search model not found! %s" % self.pdbout)
            self.error = True

    def _merge_models(self):
        """Merge all the modified models indicated at \
        :py:attr:`~swamp.searchmodel.prepare.modified_model_list` into \
        :py:attr:`~swamp.searchmodel.prepare.pdbout`"""

        if len(self.model_list) > 1:
            gesamt = Gesamt(mode='alignment', pdbin=self.modified_model_list, pdbout=self.pdbout, workdir=self.workdir,
                            logger=self.logger)
            gesamt.run()

        else:
            shutil.copyfile(self.modified_model_list[0], self.pdbout)

    def _make_workdir(self):
        """Method to crete the :py:attr:`~swamp.searchmodel.prepare.workdir` directory"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    # ------------------ Static methods ------------------

    @staticmethod
    def split_models(pdbin, directory, strip_hetatm=True):
        """Method to split an ensemble into its model components

        :param str pdbin: input pdb file name
        :param str directory: directory where the models of the ensemble will be dumped
        :param bool strip_hetatm: if set, the hetatm will be ommited from the output models (default True)
        :returns: a tuple with the output file names listed (tuple)
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

        return tuple(rslt)
