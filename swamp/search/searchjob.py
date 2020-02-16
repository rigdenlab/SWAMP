import os
import abc
import swamp
import joblib
from pyjob import Script
from swamp.wrappers import MapAlign
from swamp.logger import SwampLogger
from swamp.wrappers.aleigen import AlEigen
from swamp.utils.swamplibrary import SwampLibrary

ABC = abc.ABCMeta('ABC', (object,), {})


class SearchJob(ABC):
    """Class that implements methods to search a library of templates using contact map alignment utils
    
    The class is used to search a single subtarget against the SWAMP library. The CMO between the observed contacts
    for each member of the library and the predicted contacts of the subtarget is computed.
    
    :param id: the id given to the MrRun and identifying this job
    :type id: str, int
    :param workdir: working directory where the search will be executed and temporary files will be created
    :type workdir: str
    :param query: file name of the query contact file
    :type query: str
    :param template_library: the directory containing all the templates to be used in the search
    :type template_library: str
    :param con_format: format of the contact prediction of the query (default 'psicov')
    :type con_format: str
    :param library_format: format of the contact maps contained in the library of templates
    :type library_format: str
    :param query_pdb_benchmark: file name of the pdb structure to be used for benchmarking (default None)
    :type query_pdb_benchmark: None, str
    :param pdb_library: directory with the pdb files of the templates used for benchmarking (default None)
    :type pdb_library: str
    :param template_subset: a subset of templates to be used instead of the full library (default None)
    :type template_subset: list, tuple, None
    :param logger: logging interface to be used on the search (default None)
    :type logger: None, :obj:`swamp.logger.swamplogger.SwampLogger`
    :param algorithm: Indicate the alignment algorithm to be used in the search (default: 'mapalign')
    :type algorithm: str
    """

    def __init__(self, id, workdir, query, template_library, con_format="psicov", library_format="pdb", logger=None,
                 query_pdb_benchmark=None, pdb_library=None, template_subset=None, algorithm='mapalign',
                 python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python')):

        self._id = id
        self._query = query
        self._workdir = workdir
        self._make_workdir()
        self._template_library = template_library
        self._con_format = con_format
        self._library_format = library_format
        self._results = []
        self._query_pdb_benchmark = query_pdb_benchmark
        self._pdb_library = pdb_library
        self._template_subset = template_subset
        self._algorithm = algorithm
        self._python_interpreter = python_interpreter
        if algorithm == 'aleigen':
            self._compute_eigenvector()

        if logger is None:
            self._logger = SwampLogger(__name__)
            self.logger.init(logfile=None, use_console=True, console_level='debug')
            self.logger.info(self.search_header)
        else:
            self._logger = logger

        if self.query_pdb_benchmark is not None:
            if self.pdb_library is None:
                self.logger.warning('PDB benchmark was requested but a PDB template library was not set!')
                self.query_pdb_benchmark = None
            elif not os.path.isdir(self.pdb_library):
                self.logger.warning('PDB benchmark was requested but %s PDB library was not found!' % self.pdb_library)
                self.query_pdb_benchmark = None

    # ------------------ Some general properties ------------------

    @property
    def search_header(self):
        """Wrapper header for the logger"""

        return """**********************************************************************
********************        SWAMP - SCANJOB         ******************
**********************************************************************

"""

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, value):
        self._algorithm = value

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
    def library_format(self):
        return self._library_format

    @library_format.setter
    def library_format(self, value):
        self._library_format = value

    @property
    def python_interpreter(self):
        return self._python_interpreter

    @python_interpreter.setter
    def python_interpreter(self, value):
        self._python_interpreter = value

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
        """List of templates to be used in the search (considers the subset selection, if any)"""

        if self.template_subset is None:
            return [os.path.join(self.template_library, filename) for filename in os.listdir(self.template_library)]
        else:
            return [os.path.join(self.template_library, filename)
                    for filename in os.listdir(self.template_library)
                    if self._in_subset(filename)]

    @property
    def pickle_fname(self):
        """A pickle where the results of the search will be stored"""
        return os.path.join(self.workdir, 'search_%s_results.pckl' % self.id)

    @property
    def _python_script(self):
        """Python script to create and execute the corresponding :obj:`swamp.library.search.searchjob` instance"""

        script = 'cd {}\n{} << EOF\nfrom swamp.library.search.searchjob import SearchJob\n'.format(self.workdir,
                                                                                                   self.python_interpreter)

        attributes = ['id', 'workdir', 'query', 'template_library', 'con_format', 'library_format', 'pdb_library',
                      'query_pdb_benchmark', 'template_subset', 'algorithm']
        arguments = []

        for att in attributes:
            if self.__getattribute__(att) is not None:
                if not isinstance(self.__getattribute__(att), (list, tuple)):
                    arguments.append('%s="%s"' % (att, self.__getattribute__(att)))
                else:
                    arguments.append('%s=(%s)' % (att, ', '.join(['"%s"' % x for x in self.__getattribute__(att)])))

        script += 'job=SearchJob(%s)' % ', '.join(arguments)
        script += '\njob.run()\njob.store_pickle()\nEOF\n'

        return script

    @property
    def script(self):
        """Instance of :object:`pyjob.Script` that corresponds with the job to be executed"""

        script = Script(directory=self.workdir, prefix='searchjob_%s' % self.id, stem='', suffix='.sh')
        script.append(self._python_script)
        return script

    @property
    def _alignment_wrapper(self):
        """Location of the template library to be used during the CMO search

        :returns nothing
        :rtype None
        :raises ValueError if the algorithm is not recognised (valid algorithms: 'aleigen', 'mapalign')
        """

        if self.algorithm == 'mapalign':
            return MapAlign
        elif self.algorithm == 'aleigen':
            return AlEigen
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.algorithm)

    @property
    def eigen_template(self):
        """Eigen vector file name template"""
        return os.path.join(swamp.FRAG_EIGEN_DB, "{}.eigenvct")

    @property
    def tmp_eigen_query(self):
        """A temporary file name to contain the eigen vectors of the query"""
        return os.path.join(self.workdir, 'tmp_query.eigenvct')

    @property
    def joblist(self):
        """A list where each job corresponds with the alignment between the query and a member of the template list"""

        joblist = []

        for template in self.template_list:
            job = self._alignment_wrapper(**self._get_alignment_info(template))
            joblist.append(job)

        return joblist

    # ------------------ Protected methods ------------------

    def _get_alignment_info(self, template):
        """Create a dictionary with the arguments necessary for the contact map alignment

        :param template: the template contact map file name
        :type template: str
        :returns a dictionary with the arguments to execute a given alignment
        :rtype dict
        :raises ValueError if the algorithm is not recognised (valid algorithms: 'aleigen', 'mapalign')
        """

        if self.algorithm == 'mapalign':
            return {
                "pdb_a": self.query_pdb_benchmark,
                "pdb_b": os.path.join(self.pdb_library, self._pdbfile_template(template)),
                'workdir': self.workdir,
                'map_b': template,
                'map_a': self.query,
                'format_a': self.con_format,
                'format_b': self.library_format,
                'logger': self.logger
            }

        elif self.algorithm == 'aleigen':
            return {
                "pdb_a": self.query_pdb_benchmark,
                "pdb_b": os.path.join(self.pdb_library, self._pdbfile_template(template)),
                "workdir": self._job_dir_template(template),
                "map_b": template,
                "map_a": self.query,
                "format_b": self.library_format,
                "format_a": self.con_format,
                'logger': self.logger,
                "eigenvectors": (
                    self.tmp_eigen_query, self.eigen_template.format(os.path.basename(template).split(".")[0]))
            }

        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.algorithm)

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

    def _compute_eigenvector(self):
        """Method to create eigen vectors for the input file"""
        AlEigen.create_eigen_vectors(self.query, map_format=self.con_format, vector_output=self.tmp_eigen_query)

    # ------------------ Some general methods ------------------

    def store_pickle(self):
        """Store the list of obtained results into a pickle file"""
        joblib.dump(self.results, self.pickle_fname)

    def run(self):
        """Use the query to search against the templates using the indicated algorithm"""

        self.logger.info("Searching template library")

        for job in self.joblist:
            self.logger.debug("Alignment of query %s with template %s" % (job.map_a, job.map_b))
            job.run()
            self.results.append(job.summary_results)
            self.logger.debug(job.summary_results)

        self.logger.debug("All alignments are done!")

        if os.path.isfile(self.tmp_eigen_query):
            os.remove(self.tmp_eigen_query)
