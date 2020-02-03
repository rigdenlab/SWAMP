import abc
import os
import logging
import gemmi
from swamp.wrappers.gesamt import Gesamt

ABC = abc.ABCMeta('ABC', (object,), {})


class PrepareSearchModel(ABC):
    """Abstract class for serch model preparation"""

    def __init__(self, workdir, pdbin, pdbout, target_fa=None, logger=None):
        self._workdir = workdir
        self._pdbin = pdbin
        self._pdbout = pdbout
        self._error = False
        self._logcontents = None
        self._modified_model_list = []
        self._keywords = None
        self._target_fa = target_fa
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
    def cmd(self):
        """Abstract property to store comand to run in the terminal (if any)"""
        pass

    @property
    @abc.abstractmethod
    def modification(self):
        """Abstract property to store modification to be applied"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def target_fa(self):
        """Property to store target_fa"""
        return self._target_fa

    @target_fa.setter
    def target_fa(self, value):
        self._target_fa = value

    @property
    def logger(self):
        """Property to store logger"""
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def logcontents(self):
        """Property to store logcontents"""
        return self._logcontents

    @logcontents.setter
    def logcontents(self, value):
        self._logcontents = value

    @property
    def keywords(self):
        """Property to store keywords"""
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        self._keywords = value

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

    # ------------------ Some general methods ------------------

    def check_input(self):

        """Method to check if the required input is given"""

        if self.modification == "molrep" and self.target_fa is None:
            self.logger.error("Need to provide a target fasta sequence to use molrep!")
            self.error = True

    def check_output(self):

        """Method to check if the prepared model extists"""

        if not os.path.isfile(self.pdbout):
            self.logger.error("Modified search model not found! %s" % self.pdbout)
            self.error = True

    def _merge_models(self):

        """Create a function to merge all the modified models together"""

        gesamt = Gesamt(mode='alignment', pdbin=self.modified_model_list, pdbout=self.pdbout, workdir=self.workdir,
                        logger=self.logger)
        gesamt.run()

    def make_workdir(self):

        """Method to crete the workdir for the wrapper"""

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

    @staticmethod
    def split_models(pdbin, directory, strip_hetatm=True):

        """Method to split a pdb into its models"""

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
