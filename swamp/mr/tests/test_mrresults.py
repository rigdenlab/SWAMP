import os
import dill
import unittest
import collections
from swamp.utils import remove
from swamp.mr.mrresults import MrResults

RESULTS = collections.namedtuple('results', ['results'])
WORKDIR = os.path.join(os.environ['CCP4_SCR'], 'test_workdir')
MR_DIR = os.path.join(WORKDIR, 'swamp_mr')


class MrResultsTestCase(unittest.TestCase):

    def test_1(self):
        search_1 = os.path.join(MR_DIR, 'search_1')
        search_1_run_1 = os.path.join(MR_DIR, 'search_1', 'run_1')
        search_2 = os.path.join(MR_DIR, 'search_2')
        search_2_run_1 = os.path.join(MR_DIR, 'search_2', 'run_1')
        search_2_run_2 = os.path.join(MR_DIR, 'search_2', 'run_2')
        directories = [WORKDIR, MR_DIR, search_1, search_1_run_1, search_2, search_2_run_1, search_2_run_2]
        for directory in directories:
            if not os.path.isdir(directory):
                os.mkdir(directory)
        self.addCleanup(remove, WORKDIR)

        results = RESULTS(
            results=[['SEARCH_1', 'RUN_1', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor', 'local_CC',
                      'overall_CC', '15', 'acl', 'is_extended', 'solution']])
        with open(os.path.join(search_1_run_1, 'results.pckl'), 'wb') as fhandle:
            dill.dump(results, fhandle)
        results = RESULTS(
            results=[['SEARCH_2', 'RUN_1', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor', 'local_CC',
                      'overall_CC', '45', 'acl', 'is_extended', 'solution']])
        with open(os.path.join(search_2_run_1, 'results.pckl'), 'wb') as fhandle:
            dill.dump(results, fhandle)
        results = RESULTS(
            results=[['SEARCH_2', 'RUN_2', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor', 'local_CC',
                      'overall_CC', '9', 'acl', 'is_extended', 'solution']])
        with open(os.path.join(search_2_run_2, 'results.pckl'), 'wb') as fhandle:
            dill.dump(results, fhandle)

        results = MrResults(swamp_workdir=WORKDIR)
        self.assertTupleEqual((os.path.join(search_1_run_1, 'results.pckl'),
                               os.path.join(search_2_run_2, 'results.pckl'),
                               os.path.join(search_2_run_1, 'results.pckl')), results.pickle_list)
        self.assertListEqual([['SEARCH_2', 'RUN_2', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor',
                               'local_CC', 'overall_CC', '9', 'acl', 'is_extended', 'solution'],
                              ['SEARCH_2', 'RUN_1', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor',
                               'local_CC', 'overall_CC', '45', 'acl', 'is_extended', 'solution'],
                              ['SEARCH_1', 'RUN_1', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor',
                               'local_CC', 'overall_CC', '15', 'acl', 'is_extended', 'solution']]
                             , results.results)
        self.assertListEqual(["SEARCH ID", "RUN ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE",
                              "RFMC_RFACT", "RFMC_CC_LOC", "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "IS_EXTENDED",
                              "SOLUTION"], results._result_table_fields)
        self.assertEqual(results.logger_header, """\n******************************************************************\
****
*******************          SWAMP-MR RESULTS          ***************
**********************************************************************

Recovering results now...
""")
