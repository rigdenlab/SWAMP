import unittest
import numpy as np
from swamp.wrappers.mapalign import MapAlign


class MockMapAlign(MapAlign):
    """A class to mock :py:obj:`~swmap.wrappers.mapalign.MapAlign` for testing purposes"""

    @property
    def source(self):
        """Override :py:attr:`~swmap.wrappers.mapalign.MapAlign.source` with dummy variable"""
        return '/empty/path/map_align'

    def run(self):
        """Override :py:func:`~swmap.wrappers.mapalign.MapAlign.run` so that it does not execute \
        :py:attr:`~swmap.wrappers.mapalign.MapAlign.cmd`"""

        self.logcontents = """OPT -------------------------------------------------------------------
OPT                           MAP_ALIGN                                
OPT -------------------------------------------------------------------
OPT   -a          /home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map
OPT   -b          /home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map
OPT   -gap_o      -1
OPT   -gap_e      -0.01
OPT   -sep_cut    3
OPT   -iter       20
OPT   -silent     0
OPT -------------------------------------------------------------------
OPT   -use_gap_ss  0
OPT   -use_prf     0
OPT -------------------------------------------------------------------
TMP	0_1_0	21.7832	-1.06	20.7233
TMP	0_1_1	21.7832	-1.06	20.7233
TMP	0_1_2	24.8138	-2.525	22.2888
TMP	2_32_2	31.8905	-4.39	27.5005
TMP	2_32_3	31.8905	-4.39	27.5005
MAX 2_2_2	/home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map	/home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map	44.5272	-4.38	40.1472	74	1:1	2:2	3:3	4:4	5:5	6:6	7:7	8:8	9:9	10:10	11:11	12:12	13:13	14:14	15:20	16:21	17:22	18:23	19:2420:25	21:26	22:27	23:30	24:31	25:32	26:33	27:34	28:35	29:36	30:37	31:42	32:43	33:44	34:45	35:46	36:47	37:48	38:49	39:52	40:53	42:54	43:55	44:56	45:57	46:58	47:60	48:61	49:62	50:63	51:64	52:65	53:113	54:114	55:115	56:116	57:117	58:118	59:119	60:12061:143	62:144	63:145	64:146	65:147	66:148	67:149	68:150	69:151	70:152	71:153	72:154	73:155	74:156	75:157
"""
        self.get_scores()


class MapAlignWrapperTestCase(unittest.TestCase):

    def test_1(self):
        mapalign = MockMapAlign(workdir='/empty/path/workdir', map_a='/empty/path/map_a.psicov', format_b='psicov',
                                map_b='/empty/path/map_b.psicov', pdb_b='/empty/path/pdb_b', gap_o=-8, sep_cut=6,
                                format_a='psicov', pdb_a='/empty/path/pdb_a', gap_e=-4.6)
        self.assertListEqual(mapalign.cmd, ['/empty/path/map_align', '-a', '/empty/path/workdir/map_a.mapalign', '-b',
                                            '/empty/path/workdir/map_b.mapalign', "-sep_cut", '6', "-gap_e", '-4.6',
                                            "-gap_o", '-8', '-silent'])
        alignment = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14,
                     144: 62, 145: 63, 146: 64, 147: 65, 20: 15, 21: 16, 22: 17, 23: 18, 152: 70, 153: 71, 26: 21,
                     27: 22, 156: 74, 12061: 60, 30: 23, 31: 24, 32: 25, 33: 26, 34: 27, 35: 28, 36: 29, 37: 30,
                     116: 56, 42: 31, 43: 32, 44: 33, 45: 34, 46: 35, 47: 36, 48: 37, 49: 38, 52: 39, 53: 40, 54: 42,
                     55: 43, 56: 44, 57: 45, 58: 46, 60: 47, 61: 48, 62: 49, 63: 50, 64: 51, 65: 52, 150: 68, 155: 73,
                     151: 69, 154: 72, 157: 75, 113: 53, 114: 54, 115: 55, 2420: 19, 117: 57, 118: 58, 119: 59, 148: 66,
                     149: 67}
        mapalign.run()
        self.assertDictEqual(alignment, mapalign.alignment)
        self.assertListEqual(['/empty/path/map_a.psicov', '/empty/path/map_b.psicov', 44.5272, -4.38, 40.1472,
                              74, np.nan, np.nan, np.nan, np.nan], mapalign.summary_results)


if __name__ == '__main__':
    unittest.main()
