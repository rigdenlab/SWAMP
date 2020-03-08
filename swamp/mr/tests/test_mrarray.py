import os
import dill
import unittest
import collections
from swamp.utils import remove
from swamp.mr.mrarray import MrArray
from swamp.mr.mrjob import MrJob

RESULTS = collections.namedtuple('results', ['results'])
WORKDIR = os.path.join(os.environ['CCP4_SCR'], 'test_workdir')


class MrArrayTestCase(unittest.TestCase):

    def test_1(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440)

        self.assertTrue(os.path.isdir(WORKDIR))
        self.addCleanup(remove, WORKDIR)
        self.assertListEqual(array.cleanup_dir_list, [WORKDIR])
        self.assertDictEqual(array._other_task_info, {'directory': WORKDIR, 'shell': "/bin/bash", 'runtime': 1440,
                                                      'queue': 'dummy_queue', 'max_array_size': 100,
                                                      'environment': 'dummy_env'})
        self.assertListEqual(array._result_table_fields,
                             ["SEARCH ID", "RUN ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE",
                              "RFMC_RFACT", "RFMC_CC_LOC", "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "IS_EXTENDED",
                              "SOLUTION"])
        self.assertEqual(os.path.join(WORKDIR, "results.pckl"), array.result_pickle_fname)
        self.assertEqual(os.path.join(WORKDIR, "results.table"), array.result_table_fname)

    def test_2(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local')

        self.addCleanup(remove, WORKDIR)
        self.assertDictEqual(array._other_task_info, {'directory': WORKDIR, 'shell': "/bin/bash", 'runtime': 1440,
                                                      'queue': 'dummy_queue', 'processes': 100,
                                                      'environment': 'dummy_env'})

    def test_3(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local', phased_mtz='/empty/path/phases.mtz')
        self.addCleanup(remove, WORKDIR)

        job_1 = MrJob(id='job_1', workdir=WORKDIR, python_interpreter='/empty/path/python')
        array.add(job_1)
        job_2 = MrJob(id='job_2', workdir=WORKDIR, python_interpreter='/empty/path/python')
        array.add(job_2)
        job_3 = MrJob(id='job_3', workdir=WORKDIR, python_interpreter='/empty/path/python')
        array.add(job_3)

        self.assertEqual(job_1.phased_mtz, '/empty/path/phases.mtz')
        self.assertEqual(job_1.target_mtz, '/empty/path/target.mtz')
        self.assertEqual(job_1.target_fa, '/empty/path/target.fasta')

        job_b = MrJob(id='job_1', workdir=WORKDIR, python_interpreter='/empty/path/python')
        with self.assertRaises(ValueError):
            array.add(job_b)
        with self.assertRaises(TypeError):
            array.add('dummy-job')

        self.assertEqual(job_1, array[0])
        with self.assertRaises(NotImplementedError):
            self.assertEqual(job_1, array[:1])
        for idx, job in enumerate(array):
            self.assertEqual(job, array[idx])
        self.assertTrue(job_1.id in array)
        del array[job_2.id]
        self.assertFalse(job_2.id in array)
        self.assertEqual(2, len(array))
        self.assertListEqual([job_3, job_1], list(reversed(array)))

    def test_4(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local', phased_mtz='/empty/path/phases.mtz')

        self.addCleanup(remove, WORKDIR)

        job_1 = MrJob(id='job_1', workdir=os.path.join(WORKDIR, 'job_1'), python_interpreter='/empty/path/python')
        pickle_1_fname = os.path.join(WORKDIR, 'job_1', "results.pckl")
        results_1 = RESULTS(
            results=[['SEARCH_1', 'RUN_1', 'LLG_1', 'TFZ_1', 'local_CC_1', 'overall_CC_1', 'rfree_1', 'rfactor_1',
                      'local_CC_1', 'overall_CC_1', 'cc_1', 'acl_1', 'is_extended_1', 'solution_1']])
        with open(pickle_1_fname, 'wb') as fhandle:
            dill.dump(results_1, fhandle)
        self.addCleanup(remove, os.path.join(WORKDIR, 'job_1'))
        self.addCleanup(remove, pickle_1_fname)
        array.add(job_1)

        job_2 = MrJob(id='job_2', workdir=os.path.join(WORKDIR, 'job_2'), python_interpreter='/empty/path/python')
        self.addCleanup(remove, os.path.join(WORKDIR, 'job_2'))
        pickle_2_fname = os.path.join(WORKDIR, 'job_2', "results.pckl")
        results_2 = RESULTS(
            results=[['SEARCH_2', 'RUN_2', 'LLG_2', 'TFZ_2', 'local_CC_2', 'overall_CC_2', 'rfree_2', 'rfactor_2',
                      'local_CC_2', 'overall_CC_2', 'cc_2', 'acl_2', 'is_extended_2', 'solution_2']])
        with open(pickle_2_fname, 'wb') as fhandle:
            dill.dump(results_2, fhandle)
        self.addCleanup(remove, os.path.join(WORKDIR, 'job_2'))
        self.addCleanup(remove, pickle_2_fname)
        array.add(job_2)

        job_3 = MrJob(id='job_3', workdir=os.path.join(WORKDIR, 'job_3'), python_interpreter='/empty/path/python')
        self.addCleanup(remove, os.path.join(WORKDIR, 'job_3'))
        pickle_3_fname = os.path.join(WORKDIR, 'job_3', "results.pckl")
        results_3 = RESULTS(
            results=[['SEARCH_3', 'RUN_3', 'LLG_3', 'TFZ_3', 'local_CC_3', 'overall_CC_3', 'rfree_3', 'rfactor_3',
                      'local_CC_3', 'overall_CC_3', 'cc_3', 'acl_3', 'is_extended_3', 'solution_3']])
        with open(pickle_3_fname, 'wb') as fhandle:
            dill.dump(results_3, fhandle)
        self.addCleanup(remove, os.path.join(WORKDIR, 'job_3'))
        self.addCleanup(remove, pickle_3_fname)
        array.add(job_3)

        array.append_results()
        self.assertListEqual(array.results, [
            ['SEARCH_1', 'RUN_1', 'LLG_1', 'TFZ_1', 'local_CC_1', 'overall_CC_1', 'rfree_1', 'rfactor_1', 'local_CC_1',
             'overall_CC_1', 'cc_1', 'acl_1', 'is_extended_1', 'solution_1'],
            ['SEARCH_2', 'RUN_2', 'LLG_2', 'TFZ_2', 'local_CC_2', 'overall_CC_2', 'rfree_2', 'rfactor_2', 'local_CC_2',
             'overall_CC_2', 'cc_2', 'acl_2', 'is_extended_2', 'solution_2'],
            ['SEARCH_3', 'RUN_3', 'LLG_3', 'TFZ_3', 'local_CC_3', 'overall_CC_3', 'rfree_3', 'rfactor_3', 'local_CC_3',
             'overall_CC_3', 'cc_3', 'acl_3', 'is_extended_3', 'solution_3']])

        array.create_result_table_outfile()
        self.assertTrue(os.path.isfile(array.result_table_fname))
        self.addCleanup(remove, array.result_table_fname)
        expected_table = """+-----------+--------+-------+-------+-------------+--------------+------------+------------+-------------+--------------+---------+----------+---------------+------------+
| SEARCH ID | RUN ID |  LLG  |  TFZ  | PHSR_CC_LOC | PHSR_CC_ALL  | RFMC_RFREE | RFMC_RFACT | RFMC_CC_LOC | RFMC_CC_ALL  | SHXE_CC | SHXE_ACL |  IS_EXTENDED  |  SOLUTION  |
+-----------+--------+-------+-------+-------------+--------------+------------+------------+-------------+--------------+---------+----------+---------------+------------+
|  SEARCH_3 | RUN_3  | LLG_3 | TFZ_3 |  local_CC_3 | overall_CC_3 |  rfree_3   | rfactor_3  |  local_CC_3 | overall_CC_3 |   cc_3  |  acl_3   | is_extended_3 | solution_3 |
|  SEARCH_2 | RUN_2  | LLG_2 | TFZ_2 |  local_CC_2 | overall_CC_2 |  rfree_2   | rfactor_2  |  local_CC_2 | overall_CC_2 |   cc_2  |  acl_2   | is_extended_2 | solution_2 |
|  SEARCH_1 | RUN_1  | LLG_1 | TFZ_1 |  local_CC_1 | overall_CC_1 |  rfree_1   | rfactor_1  |  local_CC_1 | overall_CC_1 |   cc_1  |  acl_1   | is_extended_1 | solution_1 |
+-----------+--------+-------+-------+-------------+--------------+------------+------------+-------------+--------------+---------+----------+---------------+------------+"""
        with open(array.result_table_fname, 'r') as fhandle:
            actual_table = "".join(fhandle.readlines())
        self.assertEqual(expected_table, actual_table)

    def test_5(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local', phased_mtz='/empty/path/phases.mtz')

        self.addCleanup(remove, WORKDIR)
        array.logger = None
        array.store_pickle()
        self.assertTrue(os.path.isfile(array.result_pickle_fname))
        self.addCleanup(remove, array.result_pickle_fname)
        with open(array.result_pickle_fname, 'rb') as fhandle:
            test_array = dill.load(fhandle)
        self.assertEqual(vars(test_array), vars(array))

    def test_6(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local', phased_mtz='/empty/path/phases.mtz')

        self.addCleanup(remove, WORKDIR)
        self.assertTrue(os.path.isdir(WORKDIR))
        array._cleanup_files()
        self.assertFalse(os.path.isdir(WORKDIR))

    def test_7(self):
        array = MrArray(id='test', workdir=WORKDIR, target_mtz='/empty/path/target.mtz', queue_environment='dummy_env',
                        target_fa='/empty/path/target.fasta', max_concurrent_nprocs=100, queue_name='dummy_queue',
                        job_kill_time=1440, platform='local', phased_mtz='/empty/path/phases.mtz')
        self.assertListEqual(array._result_table_fields,
                             ["SEARCH ID", "RUN ID", "LLG", "TFZ", "PHSR_CC_LOC", "PHSR_CC_ALL", "RFMC_RFREE",
                              "RFMC_RFACT", "RFMC_CC_LOC", "RFMC_CC_ALL", "SHXE_CC", "SHXE_ACL", "IS_EXTENDED",
                              "SOLUTION"])
        self.assertDictEqual(array.init_params, {'id': 'test', 'job_kill_time': 1440, 'logger': None,
                                                 'max_array_size': None, 'max_concurrent_nprocs': 100,
                                                 'phased_mtz': '/empty/path/phases.mtz', 'platform': 'local',
                                                 'queue_environment': 'dummy_env', 'queue_name': 'dummy_queue',
                                                 'silent': False, 'target_fa': '/empty/path/target.fasta',
                                                 'target_mtz': '/empty/path/target.mtz',
                                                 'workdir': '/tmp/filo/test_workdir'})
        self.assertEqual(array.pipeline_header, """\n**********************************************************************
**********************           {}           ******************
**********************************************************************

""")

        self.assertEqual("""Arguments provided:

	id: test
	workdir: /tmp/filo/test_workdir
	target_mtz: /empty/path/target.mtz
	target_fa: /empty/path/target.fasta
	platform: local
	queue_name: dummy_queue
	logger: None
	max_array_size: None
	queue_environment: dummy_env
	phased_mtz: /empty/path/phases.mtz
	max_concurrent_nprocs: 100
	job_kill_time: 1440
	silent: False
""", array._inform_args(**array.init_params))