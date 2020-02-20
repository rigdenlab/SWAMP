import os
import unittest
from swamp.utils import create_tempfile
from swamp.parsers.phaserparser import PhaserParser


class PhaserParserTestCase(unittest.TestCase):

    def test_1(self):
        file_contents = """REMARK TITLE [no title set]
REMARK ENSEMBLE PDB_1 EULER  49.90  83.52 180.06 FRAC -0.288 -0.526 -0.158
CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2      8
SCALE1      0.013637 -0.000000 -0.000000        0.00000
SCALE2      0.000000  0.013637 -0.000000        0.00000
SCALE3      0.000000  0.000000  0.006115        0.00000
ATOM      1  N   ALA A   1      14.378 -30.428 -13.824  1.00 61.34           N
ATOM    210  CB  ALA A  42       2.780 -38.569 -12.366  1.00 77.20           C
END
"""

        stdout_contents = """******************************************************************************************
*** Phaser Module: PREPROCESSOR                                                  2.8.3 ***
******************************************************************************************

COMPOSITION BY ASU
COMPOSITION PROTEIN MW 34738.3 NUM 1
ELLG HIRES ON
ENSEMBLE PDB_1 PDB &
"/data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/searchmodels/searchmodel_1_polyala.pdb" RMS &
0.1 
ENSEMBLE PDB_1 DISABLE CHECK ON
HKLIN "/data1/filo/results/SWAMP_benchmarking/3zux.mtz"
JOBS 1
KEYWORDS ON
LABIN F = FP SIGF = SIGFP
MUTE ON
ROOT "3zux_phaser"
SEARCH ENSEMBLE PDB_1 
SGALT BASE  P 4w 2c
SGALT SELECT LIST
SGALT TEST P 41 2 2
XYZOUT ON
KILL TIME 1440

CPU Time: 0 days 0 hrs 0 mins 0.76 secs (      0.76 secs)
Finished: Mon Feb 17 02:38:52 2020
"""
        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = PhaserParser(stdout=stdout_contents, fname=fname)
        parser.parse()
        self.assertEqual('NA', parser.LLG)
        self.assertTrue(parser.error)

    def test_2(self):
        file_contents = """REMARK TITLE [no title set]
REMARK Log-Likelihood Gain:   70.197
REMARK  RFZ=3.0 TFZ=5.6 PAK=0 LLG=70 TFZ==6.0 LLG=70 TFZ==5.9 PAK=0 LLG=70 RFZ==3.5 TFZ==5.9
REMARK ENSEMBLE PDB_1 EULER  49.90  83.52 180.06 FRAC -0.288 -0.526 -0.158
CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2      8
SCALE1      0.013637 -0.000000 -0.000000        0.00000
SCALE2      0.000000  0.013637 -0.000000        0.00000
SCALE3      0.000000  0.000000  0.006115        0.00000
ATOM      1  N   ALA A   1      14.378 -30.428 -13.824  1.00 61.34           N
ATOM    210  CB  ALA A  42       2.780 -38.569 -12.366  1.00 77.20           C
END
"""
        stdout_contents = """******************************************************************************************
*** Phaser Module: PREPROCESSOR                                                  2.8.3 ***
******************************************************************************************

COMPOSITION BY ASU
COMPOSITION PROTEIN MW 34738.3 NUM 1
ELLG HIRES ON
ENSEMBLE PDB_1 PDB &
"/data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/searchmodels/searchmodel_1_polyala.pdb" RMS &
0.1 
ENSEMBLE PDB_1 DISABLE CHECK ON
HKLIN "/data1/filo/results/SWAMP_benchmarking/3zux.mtz"
JOBS 1
KEYWORDS ON
LABIN F = FP SIGF = SIGFP
MUTE ON
ROOT "3zux_phaser"
SEARCH ENSEMBLE PDB_1 
SGALT BASE  P 4w 2c
SGALT SELECT LIST
SGALT TEST P 41 2 2
XYZOUT ON
KILL TIME 1440

CPU Time: 0 days 0 hrs 0 mins 0.76 secs (      0.76 secs)
Finished: Mon Feb 17 02:38:52 2020



--------------
MONOMERIC ELLG
--------------

   Expected LLG (eLLG)
   -------------------
   eLLG: eLLG of ensemble alone
       eLLG   RMSD frac-scat  Ensemble
    5.29275  1.077   0.08194  PDB_1

   Resolution for eLLG target
   --------------------------
   eLLG-reso: Resolution to achieve target eLLG (225)
     eLLG-reso  Ensemble
   > 2.20(all)  PDB_1

   Resolution for eLLG target: data collection
   -------------------------------------------
   eLLG-reso: Resolution to achieve target eLLG (225) with perfect data
     eLLG-reso  Ensemble
        >1.2A   PDB_1

** Solution #1 written to MTZ file:  3zux_phaser.1.mtz
   Solution #1 annotation (history):
   SOLU SET  RFZ=3.0 TFZ=5.6 PAK=0 LLG=70 TFZ==6.0 LLG=70 TFZ==5.9 PAK=0 LLG=70 TFZ==5.9
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER   49.9   83.5  180.1 FRAC -0.29 -0.53 -0.16 BFAC -1.68 #TFZ==5.9
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48

   Annotation shown below for a maximum of a further 9 solutions
   See SOL file for any additional solutions

   Solution #2 annotation (history):
   SOLU SET  RFZ=2.8 TFZ=5.3 PAK=0 LLG=67 LLG=67 PAK=0 LLG=67
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER  241.6   60.7  143.4 FRAC -0.06 -0.40 -0.31 BFAC -1.80
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48

   Solution #3 annotation (history):
   SOLU SET  RFZ=3.0 TFZ=5.5 PAK=0 LLG=64 LLG=65 PAK=0 LLG=65
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER   49.6   79.0  184.1 FRAC -0.19  0.01 -0.16 BFAC  0.11
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48

   Solution #4 annotation (history):
   SOLU SET  RFZ=2.8 TFZ=5.7 PAK=0 LLG=63 LLG=63 PAK=0 LLG=63
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER   64.6   62.5  142.9 FRAC  0.05  0.38 -0.24 BFAC  0.22
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48

   Solution #5 annotation (history):
   SOLU SET  RFZ=3.0 TFZ=4.6 PAK=0 LLG=62 LLG=63 PAK=0 LLG=63
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER   49.9   87.8  177.8 FRAC -0.18  0.07 -0.14 BFAC  1.12
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48

   Solution #6 annotation (history):
   SOLU SET  RFZ=3.0 TFZ=4.9 PAK=8 LLG=63 LLG=63 PAK=8 LLG=63
   SOLU SPAC P 41 2 2
   SOLU 6DIM ENSE PDB_1 EULER  316.7   83.0  181.1 FRAC -0.06  0.07 -0.13 BFAC -0.01
   SOLU ENSEMBLE PDB_1 VRMS DELTA -0.3997 #RMSD  1.08  1.08  0.95  0.79 #VRMS  0.87  0.87  0.71  0.48
"""

        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = PhaserParser(stdout=stdout_contents, fname=fname)
        parser.parse()

        self.assertEqual('70.197', parser.LLG)
        self.assertEqual('5.29275', parser.eLLG)
        self.assertEqual('5.9', parser.TFZ)
        self.assertEqual('3.5', parser.RFZ)
        self.assertEqual('-0.3997', parser.VRMS)

        self.assertTupleEqual(('70.197', '5.9', '3.5', '5.29275', '-0.3997'), parser.summary)


if __name__ == '__main__':
    unittest.main()
