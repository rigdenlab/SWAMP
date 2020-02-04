import os
import swamp
from swamp.wrappers.aleigen import AlEigen
from swamp.library.scan.scan import Scan


class AlEigenScan(Scan):
    """Class that implements methods to scan a library of templates using AL-Eigen and data structures
    to store the results.

    :param kwargs: arguments when the class is instantiated are passed to swamp.scan.scan
    :type kwargs: dict
    """

    def __init__(self, **kwargs):
        super(AlEigenScan, self).__init__(**kwargs)
        self._tmp_eigen_a = AlEigen.get_tempfile()
        self.compute_eigenvector()

    @property
    def _algorithm_name(self):
        return 'aleigen'

    @property
    def _align_engine(self):
        return AlEigen

    @property
    def eigen_b_template(self):
        return os.path.join(swamp.FRAG_EIGEN_DB, "{}.eigenvct")

    @property
    def tmp_eigen_a(self):
        return self._tmp_eigen_a

    @tmp_eigen_a.setter
    def tmp_eigen_a(self, value):
        self._tmp_eigen_a = value

    @property
    def joblist(self):
        """A list of jobs to be completed. Each job corresponds with a contact map alignment"""

        joblist = []
        for template in self.template_list:
            job = {
                "pdb_a": self.query_pdb_benchmark,
                "pdb_b": os.path.join(self.pdb_library, self._pdbfile_template(template)),
                "workdir": self._job_dir_template(template),
                "map_b": template,
                "map_a": self.query,
                "format_b": self.library_format,
                "format_a": self.con_format,
                "eigenvectors": (
                    self.tmp_eigen_a, self.eigen_b_template.format(os.path.basename(template).split(".")[0]))
            }
            joblist.append(job)

        return joblist

    def compute_eigenvector(self):
        """Method to create eigen vectors for the input file"""

        AlEigen.create_eigen_vectors(self.query, map_format=self.con_format, vector_output=self.tmp_eigen_a)
