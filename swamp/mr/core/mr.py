import abc
import os
import dill
import shutil
from prettytable import PrettyTable
from swamp.logger.swamplogger import SwampLogger

ABC = abc.ABCMeta('ABC', (object,), {})


class Mr(ABC):
    """This is an abstract class for MR pipelines. It implements data structures and methods commonly used in MR tasks.

    :param str id: unique identifier for the :py:obj:`~swamp.mr.core.mr.Mr` instance
    :param str target_fa: target's fasta filename
    :param str target_mtz: target's mtz filename
    :param str workdir: working directory where the MR pipeline will be executed
    :param str phased_mtz: filename of the target's mtz containing phase information (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    :param bool silent: if set to True the logger will not print messages (default False)
    :ivar bool error: True if errors have occurred at some point on the pipeline
    :ivar list results: contains the results obtained in the MR pipeline
    """

    def __init__(self, id, target_fa, target_mtz, workdir, phased_mtz=None, logger=None, silent=False):

        self._init_params = locals()
        self._id = id
        self._search_id = id.split('_')[1]
        self._run_id = id.split('_')[3]
        self._workdir = workdir
        self.make_workdir()
        os.chdir(self.workdir)
        self._target_fa = target_fa
        self._target_mtz = target_mtz
        self._phased_mtz = phased_mtz
        self._error = False
        self._results = []
        if logger is None:
            self._logger = SwampLogger(__name__, silent=silent)
            self.logger.init(logfile=os.path.join(self.workdir, "swamp_%s.debug" % id), use_console=True,
                             console_level='info', logfile_level='debug')
        else:
            self._logger = logger
            self.logger.silent = silent

    # ------------------ Abstract methods and properties ------------------

    @abc.abstractmethod
    def run(self):
        """Abstract method to run the mr pipeline"""
        pass

    @abc.abstractmethod
    def append_results(self):
        """Abstract method to append the results obtained so far into :py:attr:`~swamp.mr.mr.Mr.results`"""
        pass

    @property
    @abc.abstractmethod
    def cleanup_dir_list(self):
        """ Abstract property to hold the location of the directories to cleanup after completion of the pipeline"""
        pass

    # ------------------ Some general properties ------------------

    @property
    def _result_table_fields(self):
        """A list of the column names in the :py:attr:`~swamp.mr.mr.Mr.table_contents`"""

        return ["SEARCH ID", "RUN ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE", "RFMC_RFACT",
                "RFMC_CC_LOC", "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "IS_EXTENDED", "SOLUTION"]

    @property
    def logger(self):
        """ :py:obj:`~swamp.logger.swamplogger.SwampLogger` logging interface for the MR pipeline (default None)"""
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def id(self):
        """Unique identifier for the MR instance. It is the result of combining :py:attr:`~swamp.mr.mr.Mr.search_id` and
         :py:attr:`~swamp.mr.mr.Mr.run_id`"""
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def search_id(self):
        """Corresponds with the unique search model identifier of the instance"""
        return self._search_id

    @search_id.setter
    def search_id(self, value):
        self._search_id = value

    @property
    def run_id(self):
        """Corresponds with the unique MR run identifier of the search model"""
        return self._run_id

    @run_id.setter
    def run_id(self, value):
        self._run_id = value

    @property
    def target_fa(self):
        """Target's fasta filename"""
        return self._target_fa

    @target_fa.setter
    def target_fa(self, value):
        self._target_fa = value

    @property
    def workdir(self):
        """working directory where the MR pipeline will be executed"""
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def phased_mtz(self):
        """Filename of the target's mtz containing phase information"""
        return self._phased_mtz

    @phased_mtz.setter
    def phased_mtz(self, value):
        """Filename of the target's mtz"""
        self._phased_mtz = value

    @property
    def target_mtz(self):
        return self._target_mtz

    @target_mtz.setter
    def target_mtz(self, value):
        self._target_mtz = value

    @property
    def error(self):
        """This property returns True if errors have been detected during the completion of the pipeline"""
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def results(self):
        """A list with the figures of merit obtained after the completion of the pipeline"""
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    @property
    def init_params(self):
        """A dictionary to store the initial parameters used to instantiate the :py:obj:`~swamp.mr.core.mr.Mr`"""

        if "self" in self._init_params:
            del self._init_params["self"]
        if '__class__' in self._init_params.keys():
            del self._init_params['__class__']
        return dict(self._init_params)

    @init_params.setter
    def init_params(self, value):
        if "self" in value:
            del value["self"]
        if "__class__" in value:
            del value['__class__']
        self._init_params = value

    @property
    def result_table_fname(self):
        """Filename where the string representation of :py:attr:`~swamp.mr.mr.Mr.table_contents` will be written"""
        return os.path.join(self.workdir, "results.table")

    @property
    def result_pickle_fname(self):
        """Filename where the :py:obj:`~swamp.mr.core.mr.Mr` instance will be pickled"""
        return os.path.join(self.workdir, "results.pckl")

    @property
    def pipeline_header(self):
        """Header displayed when initiating :py:obj:`~swamp.logger.swamplogger.SwampLogger`"""

        return """\n**********************************************************************
**********************           {}           ******************
**********************************************************************

"""

    @property
    def table_contents(self):
        """String representation of :py:attr:`~swamp.mr.mr.Mr.results` displayed as a table"""

        self.results = sorted(self.results, key=lambda x: x[10] if x[10] != 'NA' else float('0.0'), reverse=True)

        table = PrettyTable(self._result_table_fields)
        for result in self.results:
            if len(result) == len(self._result_table_fields):
                table.add_row(result)
            else:
                self.logger.warning('Results out of bounds! %s' % ', '.join(result))

        table.sortby = 'SHXE_CC'
        table.reversesort = True

        return table.get_string()

    # ------------------ Some general methods ------------------

    def make_workdir(self):
        """Method to create the working directory of the :py:obj:`~swamp.mr.core.mr.Mr` instance"""

        if not (os.path.isdir(self.workdir)):
            os.makedirs(self.workdir)

    def store_pickle(self, fname=None, mode="ab"):
        """Method to pickle the :py:obj:`~swamp.mr.core.mr.Mr` instance into a file

        :param fname: filename where the pickle will be created (default None)
        :type fname: str
        :param mode: mode that will be used to open the file handle (default 'ab')
        :type mode: str
        """

        if fname is None:
            fname = self.result_pickle_fname

        with open(fname, mode) as pickle_file:
            dill.dump(self, pickle_file, protocol=1)
        pickle_file.close()

    def create_result_table_outfile(self, fname=None):
        """Method to write string representation of :py:attr:`~swamp.mr.mr.Mr.table_contents` into
        :py:attr:`~swamp.mr.mr.Mr.result_table_fname` or into a specific file if `fname` is set.

        :param fname: filename where the result table will be created (default None)
        :type fname: str
        """

        if fname is None:
            fname = self.result_table_fname

        # If the run was aborted, nothing to log here
        if self.error:
            self.logger.warning("Previous errors prevented the creation of a result table")
            return

        # If the file doesnt exists, add the header
        if not (os.path.isfile(fname)):
            with open(fname, "w") as logfile:
                logfile.write(self.table_contents)
            logfile.close()

    # ------------------ Some protected and static methods ------------------

    def _cleanup_files(self):
        """Method to cleanup the files indicated at :py:attr:`~swamp.mr.mr.Mr.cleanup_dir_list`"""

        self.logger.info("Saving disk space now...")

        for directory in self.cleanup_dir_list:
            self.logger.debug("Removing %s" % directory)
            if os.path.isdir(directory):
                shutil.rmtree(directory)

    @staticmethod
    def _inform_args(**kwargs):
        """Create a string representation of the initial parameters used for the creation of this instance, as stored in
        :py:attr:`~swamp.mr.mr.Mr.init_params`"""

        msg = "Arguments provided:\n\n"
        for key in kwargs.keys():
            if key != "self":
                msg += "\t%s: %s\n" % (key, kwargs[key])
        return msg
