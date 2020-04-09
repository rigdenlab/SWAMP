import os
import unittest
import collections
from swamp.utils import remove
from swamp.mr.mrrun import MrRun
from mock import patch
import swamp.mr.targetdata

RESULTS = collections.namedtuple('results', ['results'])
WORKDIR = os.path.join(os.environ['CCP4_SCR'], 'test_workdir')


@patch.object(swamp.mr.targetdata.TargetData, 'get_info')
class MrRunTestCase(unittest.TestCase):

    def test_1(self, mock_my_method):
        mock_my_method.return_value = True
        mr_run = MrRun('search_1_run_2', workdir=os.path.join(os.environ['CCP4_SCR'], 'test_workdir'),
                       target_fa='/empty/path/target.fasta', target_mtz='/empty/path/target.mtz')
        self.assertTrue(os.path.isdir(os.path.join(os.environ['CCP4_SCR'], 'test_workdir')))
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'test_workdir', "searchmodels"), mr_run.searchmodel_dir)
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'test_workdir', "ideal_helices"),
                         mr_run.idealhelices_workdir)
        mr_run.target.mw = '9999'
        mr_run.target.ncopies = '1'
        mr_run.target.solvent = 0.5
        mr_run.target.nreflections = 600000
        mr_run.target.use_f = False
        mr_run.target.resolution = 2.0
        mr_run.target.seq_length = 200
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'test_workdir'))
        mr_run.logger = None
        self.assertDictEqual(
            {'early_kill': True, 'workdir': os.path.join(os.environ['CCP4_SCR'], 'test_workdir', 'phaser'),
             'timeout': 1800, 'threads': 1, 'phased_mtz': None, 'mtzfile': '/empty/path/target.mtz', 'mw': '9999',
             'packcutoff': None, 'nchains_asu': '1', 'sgalternative': 'NONE', 'peaks_rotcutoff': None, 'logger': None}
            , mr_run.phaser_info)

