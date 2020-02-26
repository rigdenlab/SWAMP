import os
import dill
import unittest
import collections
from pyjob import Script
from swamp.utils import remove
from swamp.mr.mrjob import MrJob

RESULTS = collections.namedtuple('results', ['results'])
WORKDIR = os.path.join(os.environ['CCP4_SCR'], 'test_workdir')


class MrJobTestCase(unittest.TestCase):

    def test_1(self):
        python_script = """cd %s
/empty/path/python << EOF
from swamp.mr import MrRun
mr_run = MrRun(id='test', workdir='%s', target_fa='/empty/path/target.fasta', target_mtz='/empty/path/target.mtz')
mr_run.phased_mtz = "/empty/path/phases.mtz"
mr_run.add_searchmodel(arg1="arg1", arg2="arg2", arg3="3")
if not mr_run.error:
    mr_run.run()
    mr_run.create_result_table_outfile()
    mr_run.store_pickle()
EOF
""" % (WORKDIR, WORKDIR)

        job = MrJob(id='test', workdir=WORKDIR, python_interpreter='/empty/path/python')
        self.assertTrue(os.path.isdir(WORKDIR))
        self.addCleanup(remove, WORKDIR)
        job.add_searchmodel(arg1='arg1', arg2='arg2', arg3=3)
        job.target_fa = '/empty/path/target.fasta'
        job.target_mtz = '/empty/path/target.mtz'
        job.phased_mtz = '/empty/path/phases.mtz'
        self.assertEqual(python_script, job.python_script)
        self.assertIsInstance(job.script, Script)
        self.assertEqual(python_script, job.script[0])

    def test_2(self):
        pickle_fname = os.path.join(WORKDIR, "results.pckl")
        job = MrJob(id='test', workdir=WORKDIR, python_interpreter='/empty/path/python')
        results = RESULTS(
            results=['SEARCH', 'RUN', 'LLG', 'TFZ', 'local_CC', 'overall_CC', 'rfree', 'rfactor', 'local_CC',
                     'overall_CC', 'cc', 'acl', 'is_extended', 'solution'])
        with open(pickle_fname, 'wb') as fhandle:
            dill.dump(results, fhandle)
        self.addCleanup(remove, pickle_fname)
        self.assertListEqual(results.results, job.results)
        self.addCleanup(remove, pickle_fname)

    def test_3(self):
        job = MrJob(id='test', workdir=WORKDIR, python_interpreter='/empty/path/python')
        with self.assertRaises(TypeError):
            job.parent_array = 'dummy_array'
