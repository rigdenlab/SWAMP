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
    """Class that implements methods to search a library of templates using contact map alignment tools
    
    The class is used to search a single query against the SWAMP library of templates. The CMO between the \
    observed contacts for each member of the library and the predicted contacts of the query is computed.
    
    :param str id: the id given to the MrRun and identifying this :py:obj:`~swamp.search.searchjob.SearchJob` instance
    :param str workdir: working directory where the search will be executed and temporary files will be created
    :param str query: file name of the query contact file
    :param str template_library: the directory containing all the templates to be used in the search
    :param str con_format: format of the contact prediction of the query (default 'psicov')
    :param str library_format: format of the contact maps contained in the library of templates
    :param str query_pdb_benchmark: file name of the pdb structure to be used for benchmarking (default None)
    :param str pdb_library: directory with the pdb files of the templates used for benchmarking (default None)
    :param tuple template_subset: a subset of templates to be used instead of the full library (default None)
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the search (default None)
    :param str algorithm: Indicate the alignment algorithm to be used in the search (default: 'mapalign')
    :ivar list results: a nested list with the results obtained in the search againts the library
    """

    def __init__(self, id, workdir, query, template_library, con_format="psicov", library_format="pdb", logger=None,
                 query_pdb_benchmark=None, pdb_library=None, template_subset=None, algorithm='mapalign',
                 python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python')):

        self.id = id
        self.query = query
        self.workdir = workdir
        self._make_workdir()
        self.template_library = template_library
        self.con_format = con_format
        self.library_format = library_format
        self.results = []
        self.query_pdb_benchmark = query_pdb_benchmark
        self.pdb_library = pdb_library
        self.template_subset = template_subset
        self.algorithm = algorithm
        self.python_interpreter = python_interpreter
        if algorithm == 'aleigen':
            self._compute_eigenvector()

        if logger is None:
            self.logger = SwampLogger(__name__)
            self.logger.init(logfile=None, use_console=True, console_level='debug')
            self.logger.info(self.search_header)
        else:
            self.logger = logger

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
        """Header displayed when initiating :py:obj:`~swamp.logger.swamplogger.SwampLogger`"""

        return """**********************************************************************
********************        SWAMP - SCANJOB         ******************
**********************************************************************

"""

    @property
    def template_list(self):
        """The list of templates to be used in the search (considers \
        :py:attr:`~swamp.search.searchjob.SearchJob.template_subset`)"""

        if self.template_subset is None:
            return [os.path.join(self.template_library, filename) for filename in os.listdir(self.template_library)]
        else:
            return [os.path.join(self.template_library, filename)
                    for filename in os.listdir(self.template_library)
                    if self._in_subset(filename)]

    @property
    def pickle_fname(self):
        """A pickle file where the :py:attr:`~swamp.search.searchjob.SearchJob.results` will be stored"""
        return os.path.join(self.workdir, 'search_%s_results.pckl' % self.id)

    @property
    def _python_script(self):
        """Python script to create and execute an identical :py:obj:`~swamp.search.searchjob.SearchJob` instance"""

        script = 'cd {}\n{} << EOF\nfrom swamp.search import SearchJob\n'.format(self.workdir, self.python_interpreter)

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
        """Instance of :py:obj:`pyjob.Script` that will be executed to complete the search"""

        script = Script(directory=self.workdir, prefix='searchjob_%s' % self.id, stem='', suffix='.sh')
        script.append(self._python_script)
        return script

    @property
    def _alignment_wrapper(self):
        """Pointer to the appropiate :py:obj:`~swamp.wrappers.wrapper.Wrapper` to be used in this search

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
        """A template for eigen vector file names"""
        return os.path.join(swamp.FRAG_EIGEN_DB, "{}.eigenvct")

    @property
    def tmp_eigen_query(self):
        """A temporary file name to contain the eigen vectors of :py:attr:`~swamp.search.searchjob.SearchJob.query`"""
        return os.path.join(self.workdir, 'tmp_query.eigenvct')

    @property
    def joblist(self):
        """A list of :py:obj:`~swamp.wrappers.wrapper.Wrapper` instances corresponding with the alignment between the \
        query and each member of the :py:attr:`~swamp.search.searchjob.SearchJob.template_list`"""

        joblist = []

        for template in self.template_list:
            job = self._alignment_wrapper(**self._get_alignment_info(template))
            joblist.append(job)

        return joblist

    # ------------------ Protected methods ------------------

    def _get_alignment_info(self, template):
        """Create a dictionary with the **kwargs necessary to initialise a \
        :py:obj:`~swamp.wrappers.wrapper.Wrapper` instance

        :param str template: the template contact map file name
        :returns: a dictionary with the arguments pass to :py:obj:`~swamp.wrappers.wrapper.Wrapper` (dict)
        :raises ValueError: if the algorithm is not recognised (valid algorithms: 'aleigen', 'mapalign')
        """

        result = {
            'map_b': template,
            'map_a': self.query,
            'format_a': self.con_format,
            'format_b': self.library_format,
            'logger': self.logger
        }

        if self.algorithm == 'mapalign':
            result['workdir'] = self.workdir
        elif self.algorithm == 'aleigen':
            result['workdir'] = self._job_dir_template(template)
            result["eigenvectors"] = self.tmp_eigen_query, self.eigen_template.format(os.path.basename(template).split(".")[0])
        else:
            raise ValueError("Unrecognised alignment tool! %s" % self.algorithm)

        if self.pdb_library is not None:
            result["pdb_b"] = os.path.join(self.pdb_library, self._pdbfile_template(template))
        if self.query_pdb_benchmark is not None:
            result["pdb_a"] = self.query_pdb_benchmark

        return result

    def _make_workdir(self):
        """Method to crete the :py:attr:`~swamp.search.searchjob.SearchJob.workdir`"""

        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def _pdbfile_template(self, mapfile):
        """Get the pdb file name for a given map file name"""
        if self.pdb_library is None:
            return None
        else:
            return os.path.join(self.pdb_library, "%s.pdb" % os.path.basename(mapfile).split(".")[0])

    def _job_dir_template(self, template):
        """Get the job directory for a given template name"""
        return os.path.join(self.workdir, os.path.basename(os.path.basename(template).split(".")[0]))

    def _in_subset(self, fname):
        """Check if a given fragment file name is in the :py:attr:`~swamp.search.searchjob.SearchJob.template_subset`

        :param str fname: file name of the fragment of interest
        :returns: True if the fragment is in :py:attr:`~swamp.search.searchjob.SearchJob.template_subset` (bool)
        """

        frag_id = os.path.basename(fname).split('.')[0]
        if frag_id in self.template_subset or SwampLibrary._get_reciprocal_id(frag_id) in self.template_subset:
            return True
        else:
            return False

    def _compute_eigenvector(self):
        """Create eigen vectors for :py:attr:`~swamp.search.searchjob.SearchJob.query`"""
        AlEigen.create_eigen_vectors(self.query, map_format=self.con_format, vector_output=self.tmp_eigen_query)

    # ------------------ Some general methods ------------------

    def store_pickle(self):
        """Store :py:attr:`~swamp.search.searchjob.SearchJob.results` into \
        :py:attr:`~swamp.search.searchjob.SearchJob.pickle_fname`"""
        joblib.dump(self.results, self.pickle_fname)

    def run(self):
        """Run each independent job contained at :py:attr:`~swamp.search.searchjob.SearchJob.joblist`"""

        self.logger.info("Searching template library")

        for job in self.joblist:
            self.logger.debug("Alignment of query %s with template %s" % (job.map_a, job.map_b))
            job.run()
            self.results.append(job.summary_results)
            self.logger.debug(job.summary_results)

        self.logger.debug("All alignments are done!")

        if os.path.isfile(self.tmp_eigen_query):
            os.remove(self.tmp_eigen_query)
