import os
import dill
from prettytable import PrettyTable
from argparse import ArgumentParser
from swamp.logger import SwampLogger
from swamp.command_line import check_path_exists


def create_argument_parser():
    """Create a parser for the command line arguments used in swamp-results"""

    parser = ArgumentParser()
    parser.add_argument("swamp_workdir", type=check_path_exists,
                        help="The workding directory where SWAMP results can be extracted")
    parser.add_argument("-top_results", type=int, nargs="?", default=50, help='Number of top results to display')
    parser.add_argument("-sort_by", type=str, nargs="?", default='SHXE_CC', help='Field to use for sorting results')

    return parser


class MrResults(object):
    """Class to handle the results py:obj:`~swamp.mr.mrarray.MrArray` instance

    :param str swamp_workdir: the :py:obj:`~swamp.mr.mrarray.MrArray` working directory
    :param `~swamp.logger.swamplogger.SwampLogger` logger: logging interface for the MR pipeline (default None)
    :ivar list results: nested list with the figures of merit obtained after the completion of \
    :py:obj:`~swamp.mr.mrarray.MrArray` isntance
    :ivar list pickle_list: a list of pickle file names for all the :py:obj:`~swamp.mr.mrrun.MrRun` instances \
    contained in the :py:obj:`~swamp.mr.mrarray.MrArray` instance
    :ivar str mr_dir: the directory with the MR results obtained with :py:obj:`~swamp.mr.mrarray.MrArray` instance

    """

    def __init__(self, swamp_workdir, logger=None):

        assert os.path.isdir(swamp_workdir), "Cannot find work_dir: %s" % swamp_workdir
        self.mr_dir = os.path.join(swamp_workdir, 'swamp_mr')
        assert os.path.isdir(self.mr_dir), "Cannot find Mr results directory: %s" % self.mr_dir
        self.pickle_list = self.lookup_pickles(self.mr_dir)
        if logger is None:
            self.logger = SwampLogger(__name__)
            self.logger.init(use_console=True, console_level='info')
        else:
            self.logger = logger

        self.logger.info(self.logger_header)
        self.results = []
        self._recover_results()

    @property
    def logger_header(self):
        """Header used whe initialising the :py:obj:`~swamp.logger.swamplogger.SwampLogger`"""

        return """\n**********************************************************************
*******************          SWAMP-MR RESULTS          ***************
**********************************************************************

Recovering results now...
"""

    @property
    def _result_table_fields(self):
        """List of the field names in the results table printed at \
        :py:func:`~swamp.mr.mrresults.MrResults.report_results`"""

        return ["SEARCH ID", "RUN ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE", "RFMC_RFACT",
                "RFMC_CC_LOC", "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "IS_EXTENDED", "SOLUTION"]

    def report_results(self, top=50, sort_by='SHXE_CC'):
        """Print in the stdout a table with the indicated top :py:attr:`~swamp.mr.mrresults.MrResults.results`

        :param int top: top number of results to display
        :param str sort_by: the column name to use for sorting the table
        """

        table = PrettyTable(self._result_table_fields)
        for result in self.results[:top]:
            if len(result) == len(self._result_table_fields):
                table.add_row(result)
            else:
                self.logger.warning('Results out of bounds! %s' % ', '.join(result))

        table.sortby = sort_by
        table.reversesort = True

        self.logger.info('Results retrieved. Showing top %s scoring results:\n\n%s' % (top, table.get_string()))

    def _recover_results(self):
        """Recover the results stored in each pickle at :py:attr:`~swamp.mr.mrresults.MrResults.pickle_list`"""

        for pickle_fname in self.pickle_list:
            with open(pickle_fname, 'rb') as pickle_fhandle:
                mr_run = dill.load(pickle_fhandle)
            pickle_fhandle.close()
            self.results += mr_run.results

        self.results = sorted(self.results, key=lambda x: x[10] if x[10] != 'NA' else float('0.0'), reverse=True)

    @staticmethod
    def lookup_pickles(mr_dir):
        """Search a given directory for result pickle files

        :param str mr_dir: the directory of interest to be searched
        :returns: a tuple with the pickle files that were found
        """

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
        """Check if a given directory is a swamp-mr search directory

        :param str dir: directory of interest
        :returns: True if the input directory is a swamp search directory otherwise False
        """
        if os.path.isdir(dir) and 'search_' in dir:
            return True
        else:
            return False

    @staticmethod
    def is_rundir(dir):
        """Check if a given directory is a valid swamp-mr run directory

        :param str dir: directory of interest
        :returns: True if the input directory is an existing swamp run directory otherwise False

        """
        if os.path.isdir(dir) and 'run_' in dir:
            return True
        else:
            return False


if __name__ == "__main__":
    parser = create_argument_parser()
    args = parser.parse_args()
    results = MrResults(args.swamp_workdir)
    results.report_results(args.top_results)
