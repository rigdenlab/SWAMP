import os
import unittest
import numpy as np
from swamp.wrappers.aleigen import AlEigen


class MockAleigen(AlEigen):
    """A class to mock :py:obj:`~swmap.wrappers.aleigen.AlEigen` for testing purposes"""

    def run(self):
        """Override :py:func:`~swmap.wrappers.aleigen.AlEigen.run` so that it does not execute \
        :py:attr:`~swmap.wrappers.aleigen.AlEigen.cmd`"""

        self.logcontents = """Score    C1    C2    CMO
0.073469 101   144   9    

11 0
12 1
15 2
16 3
18 4
19 5
20 6
21 7
22 8
23 9
24 10
25 11
26 12
27 13
29 14
30 15
33 16
34 17
37 18
52 19
53 20
54 22
55 23
56 24
57 25
59 26
60 27
63 28
64 29
66 30
67 31
72 32
74 33
75 34
76 35
78 36
79 37
82 38
83 39
85 40
86 41
89 42
93 43
94 178
95 179
96 180
97 181
98 182
99 183
100 184
"""
        self.get_scores()


class AlEigenWrapperTestCase(unittest.TestCase):

    def test_1(self):
        aleigen = MockAleigen(workdir='/empty/path/workdir', map_a='/empty/path/map_a.psicov', format_b='psicov',
                              map_b='/empty/path/map_b.psicov', pdb_b='/empty/path/pdb_b', format_a='psicov',
                              eigenvectors=('/empty/path/eigen_a', '/empty/path/eigen_b'), get_matrix=True,
                              n_eigen="15", pdb_a='/empty/path/pdb_a')
        self.assertListEqual(aleigen.cmd, [os.path.join(os.environ['CCP4'], 'bin', 'aleigen'),
                                           '/empty/path/workdir/map_a.aleigen', '/empty/path/workdir/map_b.aleigen',
                                           '15', '-e', '/empty/path/eigen_a', '/empty/path/eigen_b', '-m'])
        alignment = {0: 11, 1: 12, 2: 15, 3: 16, 4: 18, 5: 19, 6: 20, 7: 21, 8: 22, 9: 23, 10: 24, 11: 25, 12: 26,
                     13: 27, 14: 29, 15: 30, 16: 33, 17: 34, 18: 37, 19: 52, 20: 53, 22: 54, 23: 55, 24: 56, 25: 57,
                     26: 59, 27: 60, 28: 63, 29: 64, 30: 66, 31: 67, 32: 72, 33: 74, 34: 75, 35: 76, 36: 78, 37: 79,
                     38: 82, 39: 83, 40: 85, 41: 86, 42: 89, 43: 93, 178: 94, 179: 95, 180: 96, 181: 97, 182: 98,
                     183: 99, 184: 100}
        aleigen.error = True
        aleigen.run()
        self.assertDictEqual(alignment, aleigen.alignment)
        self.assertListEqual(['/empty/path/map_a.psicov', '/empty/path/map_b.psicov', 0.073469, 101,
                              144, 9.0, len(alignment), np.nan, np.nan, np.nan, np.nan], aleigen.summary_results)


if __name__ == '__main__':
    unittest.main()
