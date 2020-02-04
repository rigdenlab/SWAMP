import os
import abc
import threading
import logging
from swamp.library_tools.swamplibrary import SwampLibrary
from swamp.logger.swamplogger import SwampLogger
from swamp.scan.threading_results import Results

ABC = abc.ABCMeta('ABC', (object,), {})


class Scan(ABC):
    """Class that implements methods to scan a library of templates using contact map alignment tools

    :param workdir: working directory where the scan will be executed and temporary files will be created
    :type workdir: str
    :param query: file name of the query contact file
    :type query: str
    :param template_library: the directory containing all the templates to be used in the scan
    :type template_library: str
    :param con_format: format of the contact prediction of the query (default 'psicov')
    :type con_format: str
    :param nthreads: number of threads to be used in the scan (default 1)
    :type nthreads: int
    :param library_format: format of the contact maps contained in the library of templates
    :type library_format: str
    :param query_pdb_benchmark: file name of the pdb structure to be used for benchmarking (default None)
    :type query_pdb_benchmark: None, str
    :param pdb_library: directory with the pdb files of the templates used for benchmarking (default None)
    :type pdb_library: str
    :param template_subset: a subset of templates to be used instead of the full library (default None)
    :type template_subset: list, tuple, None
    :param logger: logging interface to be used on the scan (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`
    """

    def __init__(self, workdir, query, template_library, con_format="psicov", nthreads=1, library_format="pdb",
                 query_pdb_benchmark=None, pdb_library=None, template_subset=None, logger=None):
        self._query = query
        self._workdir = workdir
        self._make_workdir()
        self._template_library = template_library
        self._con_format = con_format
        self._nthreads = nthreads
        self._library_format = library_format
        self._logger = None
        self._initiate_logger()
        self._semaphore = threading.Semaphore(value=self.nthreads)
        self._results = Results(self.logger)
        self._child_threads = None
        self._query_pdb_benchmark = query_pdb_benchmark
        self._pdb_library = pdb_library
        self._template_subset = template_subset
        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(logfile=None, use_console=True, console_level='info')
            self.logger.info(self.scan_header)
        else:
            self._logger = logger

    # ------------------ Abstract properties and methods ------------------

    @property
    @abc.abstractmethod
    def _algorithm_name(self):
        """Abstract property to hold the name of the algorithm for the alignment"""
        pass

    @property
    @abc.abstractmethod
    def _align_engine(self):
        """Abstract property with the alignment algorithm to be used in the scan"""
        pass

    @property
    @abc.abstractmethod
    def joblist(self):
        """Abstract property with the list f jobs to run"""
        pass

        # ------------------ Some general properties ------------------

    @property
    def scan_header(self):
        """Wrapper header for the logger"""

        return """**********************************************************************
**********************           %s           ******************
**********************************************************************

""" % self._algorithm_name.upper()

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def template_subset(self):
        return self._template_subset

    @template_subset.setter
    def template_subset(self, value):
        self._template_subset = value

    @property
    def pdb_library(self):
        return self._pdb_library

    @pdb_library.setter
    def pdb_library(self, value):
        self._pdb_library = value

    @property
    def query_pdb_benchmark(self):
        return self._query_pdb_benchmark

    @query_pdb_benchmark.setter
    def query_pdb_benchmark(self, value):
        self._query_pdb_benchmark = value

    @property
    def child_threads(self):
        return self._child_threads

    @child_threads.setter
    def child_threads(self, value):
        self._child_threads = value

    @property
    def semaphore(self):
        return self._semaphore

    @semaphore.setter
    def semaphore(self, value):
        self._semaphore = value

    @property
    def library_format(self):
        return self._library_format

    @library_format.setter
    def library_format(self, value):
        self._library_format = value

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        self._results = value

    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = value

    @property
    def con_format(self):
        return self._con_format

    @con_format.setter
    def con_format(self, value):
        self._con_format = value

    @property
    def nthreads(self):
        return self._nthreads

    @nthreads.setter
    def nthreads(self, value):
        self._nthreads = value

    @property
    def template_library(self):
        return self._template_library

    @template_library.setter
    def template_library(self, value):
        self._template_library = value

    @property
    def query(self):
        return self._query

    @query.setter
    def query(self, value):
        self._query = value

    @property
    def template_list(self):
        """List of templates to be used in the scan (considers the subset selection, if any)"""

        if self.template_subset is None:
            return [os.path.join(self.template_library, filename) for filename in os.listdir(self.template_library)]
        else:
            return [os.path.join(self.template_library, filename)
                    for filename in os.listdir(self.template_library)
                    if self._in_subset(filename)]

    # ------------------ Protected methods ------------------

    def _make_workdir(self):
        """Method to crete the workdir for the wrapper"""

        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def _pdbfile_template(self, mapfile):
        """Get the pdb file name for a given map file"""

        return os.path.join(self.pdb_library, "%s.pdb" % os.path.basename(mapfile).split(".")[0])

    def _job_dir_template(self, template):
        """Get the job directory for a given template"""

        return os.path.join(self.workdir, os.path.basename(os.path.basename(template).split(".")[0]))

    def _initiate_logger(self):
        """Method to initiate the scan run logger"""

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def _in_subset(self, fname):
        """Check if a given fragment file is in the user specified subset

        :param fname: file name of the fragment of interest
        :type fname: str
        :returns True if the fragment is in the subset
        :rtype bool
        """

        frag_id = os.path.basename(fname).split('.')[0]
        if frag_id in self.template_subset or SwampLibrary._get_reciprocal_id(frag_id) in self.template_subset:
            return True
        else:
            return False

    def _align(self, job):
        """Compute the contact map alignment using multiple threads were the execution is controlled by a semaphore

        :param job: the job arguments to be used
        :type job: dict
        :returns nothing
        :rtype None
        """

        self.semaphore.acquire()
        self.logger.debug("Starting alignment of query %s with template %s" % (job["map_a"], job["map_b"]))
        my_run = self._align_engine(**job)
        my_run.run()
        self.results.register(my_run.summary_results)
        self.logger.debug(my_run.summary_results)
        self.logger.debug("Done with alignment of query %s with template %s" % (job["map_a"], job["map_b"]))
        self.semaphore.release()

    def _join_threads(self):
        """Join all the child threads so that the scanner waits for their completion"""

        mainthread = threading.current_thread()
        self.logger.debug("Mainthread is %s" % mainthread.getName())
        for t in threading.enumerate():
            if t is not mainthread and t.getName() in self.child_threads:
                self.logger.debug("joining thread %s" % t.getName())
                t.join()

    # ------------------ Some general methods ------------------

    def run(self):
        """Scan the library of templates using the indicated algorithm"""

        self.child_threads = []
        self.logger.info("Scanning template library using %s parallel threads" % self.nthreads)
        for job in self.joblist:
            t = threading.Thread(target=self._align, args=(job,))
            self.logger.debug("Sending thread %s" % t.getName())
            self.child_threads.append(t.getName())
            t.start()
        self.logger.info("Waiting for workers...\n")
        self._join_threads()
        self.logger.debug("Workers are done!")

        if hasattr(self, 'tmp_eigen_a'):
            os.remove(self.tmp_eigen_a)
