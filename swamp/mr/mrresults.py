import os
import dill
from prettytable import PrettyTable
from argparse import ArgumentParser
from swamp.logger.swamplogger import SwampLogger


class MrResults(object):
    """Class to handle the results from a SWAMP MR array task

    :param workdir: working directory where swamp-mr was executed
    :type workdir: str
    """

    def __init__(self, workdir):

        assert os.path.isdir(workdir), "Cannot find work_dir: %s" % workdir
        self.mr_dir = os.path.join(workdir, 'swamp_mr')
        assert os.path.isdir(self.mr_dir), "Cannot find Mr results directory: %s" % self.mr_dir
        self.pickle_list = self.lookup_pickles(self.mr_dir)
        self.logger = SwampLogger(__name__)
        self.logger.init(use_console=True, console_level='info')
        self.logger.info(self.logger_header)
        self.results = []
        self._recover_results()

    @property
    def logger_header(self):
        """Header used at the top of the logger"""

        return """\n**********************************************************************
*******************          SWAMP-MR RESULTS          ***************
**********************************************************************

Recovering results now...
"""

    @property
    def _result_table_fields(self):
        """List of the field names in the results table"""

        return ["SEARCH ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE", "RFMC_RFACT", "RFMC_CC_LOC",
                "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "ANOMALOUS", "CRNK2_RFREE", "CRNK2_RFACT", "IS_EXTENDED",
                "SOLUTION"]

    def report_results(self, top=50):
        """Print in the stdout the table with the top results

        :param top: top number of results to display
        :type top: int
        """

        table = PrettyTable(self._result_table_fields)
        for result in self.results[:top]:
            if len(result) == len(self._result_table_fields):
                table.add_row(result)
            else:
                self.logger.warning('Results out of bounds! %s' % ', '.join(result))

        table.sortby = 'SHXE_CC'
        table.reversesort = True

        self.logger.info('Results retrieved. Showing top %s scoring results:\n\n%s' % (top, table.get_string()))

    def _recover_results(self):
        """Recover the results stored in each Mr run pickle"""

        for pickle_fname in self.pickle_list:
            with open(pickle_fname, 'rb') as pickle_fhandle:
                mr_run = dill.load(pickle_fhandle)
            pickle_fhandle.close()
            self.results += mr_run.results

        self.results = sorted(self.results, key=lambda x: x[9] if x[9] != 'NA' else float('0.0'), reverse=True)

    @staticmethod
    def lookup_pickles(mr_dir):
        """Search a given mr directory for result pickle files"""

        pickle_list = []
        for search_dir in filter(MrResults.is_searchdir, [os.path.join(mr_dir, dir) for dir in os.listdir(mr_dir)]):
            search_dir = os.path.join(mr_dir, search_dir)
            for run_dir in filter(MrResults.is_rundir,
                                  [os.path.join(search_dir, dir) for dir in os.listdir(search_dir)]):
                run_dir = os.path.join(search_dir, run_dir)
                pickle_fname = os.path.join(run_dir, 'results.pckl')
                if os.path.isfile(pickle_fname):
                    pickle_list.append(pickle_fname)
        return tuple(pickle_list)

    @staticmethod
    def is_searchdir(dir):
        """Check if a given directory is a mr swamp search directory

        :param dir: directory of interest
        :type dir: str
        :returns True if the input directory is a swamp search directory
        :rtype bool

        """
        if os.path.isdir(dir) and 'search_' in dir:
            return True
        else:
            return False

    @staticmethod
    def is_rundir(dir):
        """Check if a given directory is a mr swamp run directory

        :param dir: directory of interest
        :type dir: str
        :returns True if the input directory is a swamp run directory
        :rtype bool

        """
        if os.path.isdir(dir) and 'run_' in dir:
            return True
        else:
            return False


if __name__ == "__main__":
    # User input
    parser = ArgumentParser()
    parser.add_argument(dest="swamp_workdir")
    parser.add_argument("-top_results", type=int, nargs="?", default=50, help='Number of top results to display')
    args = parser.parse_args()

    # Recover and display results
    results = MrResults(args.swamp_workdir)
    results.report_results(args.top_results)
