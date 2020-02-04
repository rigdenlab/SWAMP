import os
from swamp.wrappers.mapalign import MapAlign
from swamp.scan.scan import Scan


class MapAlignScan(Scan):
    """Class that implements methods to scan a library of templates using map_align and data structures
    to store the results.

    :example

    >>> my_scanner = MapAlignScan(workdir='<workdir>', query='<query_con>', con_format='<query_format>',
    >>>                      template_library='<library>', nthreads='<nthreads>', library_format='<library_format>',
    >>>                      gesamt_benchmark=True, pdb_query='<pdb_query>', pdb_library='<library_pdb>')
    >>> my_scanner.run()
    """

    @property
    def _algorithm_name(self):
        return 'mapalign'

    @property
    def _align_engine(self):
        return MapAlign

    @property
    def joblist(self):
        """A list of jobs to be completed. Each job corresponds with a contact map alignment"""

        joblist = []
        for template in self.template_list:
            job = {
                "pdb_a": self.query_pdb_benchmark,
                "pdb_b": os.path.join(self.pdb_library, self._pdbfile_template(template)),
                'workdir': self._job_dir_template(template),
                'map_b': template,
                'map_a': self.query,
                'format_a': self.con_format,
                'format_b': self.library_format
            }
            joblist.append(job)

        return joblist
