import os
import unittest
from swamp.utils import create_tempfile
from swamp.mr.targetdata import TargetData


class TargetDataTestCase(unittest.TestCase):

    def test_1(self):
        mw = 27480.2
        volume = 101412.70945162437
        self.assertTupleEqual(TargetData.estimate_contents(volume, mw), (1, 0.7))
        mw = 34738.3
        volume = 109911.785116
        self.assertTupleEqual(TargetData.estimate_contents(volume, mw), (1, 0.6))

        mw = 21476.76
        volume = 40984.630560366
        self.assertTupleEqual(TargetData.estimate_contents(volume, mw), (1, 0.4))

        mw = 22747.72
        volume = 101488.40414591249
        self.assertTupleEqual(TargetData.estimate_contents(volume, mw), (1, 0.7))

    def test_2(self):
        file_contents = """>4RI2:A|PDBID|CHAIN|SEQUENCE
LFKSKAKAPKKVEKPKLKVEDGLFGTSGGIGFTKENELFVGRVAMIGFAASLLGEGITGKGILSQLNLETGIPIYEAEPL
LLFFILFTLLGAIGALGDRGRFVDEPTTGLEKAVIPPGKDVRSALGLKTKGPLFGFTKSNELFVGRLAQLGFAFSLIGEI
ITGKGALAQLNIETGVPINEIEPLVLLNVVFFFIAAINPGTGKFITDDEEED
>4RI2:B|PDBID|CHAIN|SEQUENCE
LFKSKAKAPKKVEKPKLKVEDGLFGTSGGIGFTKENELFVGRVAMIGFAASLLGEGITGKGILSQLNLETGIPIYEAEPL
LLFFILFTLLGAIGALGDRGRFVDEPTTGLEKAVIPPGKDVRSALGLKTKGPLFGFTKSNELFVGRLAQLGFAFSLIGEI
ITGKGALAQLNIETGVPINEIEPLVLLNVVFFFIAAINPGTGKFITDDEEED
"""
        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        mw, seq_length = TargetData.read_fasta(fname)
        self.assertEqual(mw, 22526.01)
        self.assertEqual(seq_length, 212)
