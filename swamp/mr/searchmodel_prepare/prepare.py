import abc
import os
import gemmi
import shutil
import logging
from swamp.wrappers.gesamt import Gesamt

ABC = abc.ABCMeta('ABC', (object,), {})


class PrepareSearchModel(ABC):
    """Abstract class for serch model preparation

    :param workdir: working directory where the search model will be prepared
    :type workdir: str
    :param pdbin: input pdb file name
    :type pdbin: str
    :param pdbout: output pdb file name
    :type pdbout: str
    :param logger: logging interface for the instance (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`
    """

    def __init__(self, workdir, pdbin, pdbout, logger=None):
        self._workdir = workdir
        self._pdbin = pdbin
        self._pdbout = pdbout
        self._error = False
        self._logcontents = None
        self._modified_model_list = []
        if logger is None:
            self._logger = logging.getLogger(__name__)
        else:
            self._logger = logger

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def prepare(self):
        """ Abstract method to run the wrapper"""
        pass

    @property
    @abc.abstractmethod
    def modification(self):
        """Abstract property to store modification to be applied"""
        pass

    # ------------------ Prroperties ------------------

    @property
    def logger(self):
        """Property to store logger"""
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def modified_model_list(self):
        """Property to store modified_model_list"""
        return self._modified_model_list

    @modified_model_list.setter
    def modified_model_list(self, value):
        self._modified_model_list = value

    @property
    def model_list(self):
        """Property to store model_list"""
        if not os.path.isdir(self.model_dir):
            os.makedirs(self.model_dir)
            return self.split_models(pdbin=self.pdbin, directory=self.model_dir, strip_hetatm=True)
        else:
            return [os.path.join(self.model_dir, filename) for filename in os.listdir(self.model_dir)]

    @property
    def pdbout(self):
        """Property to store pdbout"""
        return self._pdbout

    @pdbout.setter
    def pdbout(self, value):
        self._pdbout = value

    @property
    def workdir(self):
        """Property to store workdir"""
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def pdbin(self):
        """Property to store pdbin"""
        return self._pdbin

    @pdbin.setter
    def pdbin(self, value):
        self._pdbin = value

    @property
    def error(self):
        """Property to store the presence of an error"""
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def model_dir(self):
        return os.path.join(self.workdir, "models")

    @property
    def _modified_model_template(self):
        """Propety to conain the modified model name template"""
        return os.path.join(self.workdir, "{}_%s.pdb" % self.modification)

    @property
    def _tmp_pdb(self):
        """Temporary pdb file to be created if necessary"""
        return os.path.join(self.workdir, "tmp_out_%s.pdb" % self.modification)

    # ------------------ Hidden methods ------------------

    def _check_output(self):
        """Check if the output file has been created, set :att:`error` if not"""

        if not os.path.isfile(self.pdbout):
            self.logger.error("Modified search model not found! %s" % self.pdbout)
            self.error = True

    def _merge_models(self):
        """Create a function to merge all the modified models together into an ensemble"""

        if len(self.model_list) > 1:
            gesamt = Gesamt(mode='alignment', pdbin=self.modified_model_list, pdbout=self.pdbout, workdir=self.workdir,
                            logger=self.logger)
            gesamt.run()

        else:
            shutil.copyfile(self.modified_model_list[0], self.pdbout)

    def _make_workdir(self):
        """Method to crete the workdir for the search model preparation"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    # ------------------ Static methods ------------------

    @staticmethod
    def split_models(pdbin, directory, strip_hetatm=True):
        """Method to split an ensemble into its model components

        :param pdbin: input pdb file name
        :type pdbin: str
        :param directory: directory where the models of the ensemble will be dumped
        :type directory: str
        :param strip_hetatm: if set, the hetatm will be ommited from the output models (default True)
        :type strip_hetatm: bool
        :return a tuple with the output file names listed
        :rtype tuple
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
