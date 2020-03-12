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

    def test_3(self):
        file_contents = """REMARK TITLE [no title set]
REMARK  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.1 TFZ=7.4 PAK=39 LLG=318 TFZ==10.2
REMARK ENSEMBLE FIX EULER   0.00   0.00   0.00 FRAC -0.000 -0.000 0.000
REMARK ENSEMBLE PDB_idealhelix EULER  16.22 122.55  39.80 FRAC -0.159 0.068 -0.027
REMARK ENSEMBLE PDB_idealhelix EULER  14.97  62.84 219.41 FRAC -0.247 0.202 -0.328
REMARK ENSEMBLE PDB_idealhelix EULER 139.16 107.74  52.50 FRAC 0.095 -0.076 -0.109
REMARK ENSEMBLE PDB_idealhelix EULER  72.15  80.04 213.17 FRAC -0.040 -0.258 -0.353
CRYST1   49.754   72.560   95.775  90.00  90.00  90.00 P 21 21 21   16
SCALE1      0.020099 -0.000000 -0.000000        0.00000
SCALE2      0.000000  0.013782 -0.000000        0.00000
SCALE3      0.000000  0.000000  0.010441        0.00000
ATOM      1  N   ALA A   1      -3.826  13.790 -31.009  1.00 32.20           N
"""
        stdout_contents = """"******************************************************************************************
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
       12.0  0.305   0.04195  PDB_idealhelix

   Resolution for eLLG target
   --------------------------
   eLLG-reso: Resolution to achieve target eLLG (225)
     eLLG-reso  Ensemble
   > 2.10(all)  PDB_idealhelix

   Resolution for eLLG target: data collection
   -------------------------------------------
   eLLG-reso: Resolution to achieve target eLLG (225) with perfect data
     eLLG-reso  Ensemble
        >1.2A   PDB_idealhelix

------------
OUTPUT FILES
------------

   No files output

CPU Time: 0 days 3 hrs 10 mins 25.44 secs (  11425.44 secs)
Finished: Thu Mar 12 04:02:50 2020

******************************************************************************************
*** Phaser Module: AUTOMATED MOLECULAR REPLACEMENT                               2.8.3 ***
******************************************************************************************

** No solutions after final packing
** Solution reverts to previous partial solution
**    Only packing solutions retained
**       34 solution(s) retained of 34
** Number of solutions = 34
** Solution TOP LLG = 317.6
** Solution TOP TFZ = 10.2

CPU Time: 0 days 3 hrs 10 mins 25.44 secs (  11425.44 secs)
Finished: Thu Mar 12 04:02:50 2020

******************************************************************************************
*** Phaser Module: AUTOMATED MOLECULAR REPLACEMENT                               2.8.3 ***
******************************************************************************************

** Sorry - No solution with all components
   You may find a solution with a different search or selection strategy

** Solutions written to SOL file:  5hxc_phaser.sol
** Pdb and/or Mtz files have been written with results for 1 of these solutions

** Partial Solution #1 written to PDB file:  5hxc_phaser.1.pdb
** Partial Solution #1 written to PDB (ensemble) file:  5hxc_phaser.1.1.pdb
** Partial Solution #1 written to PDB (ensemble) file:  5hxc_phaser.1.2.pdb
** Partial Solution #1 written to PDB (ensemble) file:  5hxc_phaser.1.3.pdb
** Partial Solution #1 written to PDB (ensemble) file:  5hxc_phaser.1.4.pdb
** Partial Solution #1 written to PDB (ensemble) file:  5hxc_phaser.1.5.pdb
   Partial Solution #1 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.1
    TFZ=7.4 PAK=39 LLG=318 TFZ==10.2
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   72.2   80.0  213.2 FRAC -0.04 -0.26 -0.35 BFAC  7.09 #TFZ==10.2
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Annotation shown below for a maximum of a further 9 solutions
   See SOL file for any additional solutions

   Partial Solution #2 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.2
    TFZ=8.2 PAK=39 LLG=314 TFZ==8.9
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  277.0   43.6  223.4 FRAC -0.04 -0.20 -0.34 BFAC  8.82 #TFZ==8.9
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #3 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.0
    TFZ=6.8 PAK=39 LLG=314 TFZ==8.9
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  340.6   51.2  202.6 FRAC -0.06 -0.24 -0.38 BFAC  8.91 #TFZ==8.9
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #4 annotation (history):
   SOLU SET  PAK=0 RFZ=4.3 TFZ=5.4 PAK=38 LLG=156 RFZ=3.4 TFZ=7.1 PAK=38 LLG=208 TFZ==7.3 RFZ=4.0 TFZ=8.5 PAK=38
    RFZ=3.1 TFZ=7.3 PAK=38 LLG=313 TFZ==9.9
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.3   62.2  219.9 FRAC -0.24  0.20 -0.33 BFAC  6.17
   SOLU 6DIM ENSE PDB_idealhelix EULER  275.9  114.0   38.4 FRAC -0.14  0.10 -0.05 BFAC  6.88 #TFZ==7.3
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.09 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   72.3   80.2  213.1 FRAC -0.04 -0.26 -0.35 BFAC  7.51 #TFZ==9.9
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #5 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.2
    TFZ=6.8 PAK=39 LLG=311 TFZ==9.1
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  146.3  109.4   22.0 FRAC  0.06  0.22 -0.08 BFAC 10.29 #TFZ==9.1
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #6 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=2.9
    TFZ=7.8 PAK=39 LLG=310 TFZ==8.6
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   58.1   99.6   52.4 FRAC  0.01  0.22 -0.10 BFAC  9.77 #TFZ==8.6
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #7 annotation (history):
   SOLU SET  PAK=0 RFZ=4.3 TFZ=5.4 PAK=38 LLG=156 RFZ=3.4 TFZ=7.1 PAK=38 LLG=208 TFZ==7.3 RFZ=4.0 TFZ=8.5 PAK=38
    RFZ=3.0 TFZ=6.4 PAK=38 LLG=309 TFZ==8.6
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.3   62.2  219.9 FRAC -0.24  0.20 -0.33 BFAC  6.17
   SOLU 6DIM ENSE PDB_idealhelix EULER  275.9  114.0   38.4 FRAC -0.14  0.10 -0.05 BFAC  6.88 #TFZ==7.3
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.09 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  340.5   51.2  202.5 FRAC -0.06 -0.24 -0.38 BFAC  9.53 #TFZ==8.6
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #8 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.1
    TFZ=6.3 PAK=39 LLG=309 TFZ==8.6
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   46.4   51.7  229.8 FRAC  0.17  0.14 -0.24 BFAC 13.49 #TFZ==8.6
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #9 annotation (history):
   SOLU SET  PAK=0 RFZ=4.3 TFZ=5.4 PAK=38 LLG=156 RFZ=3.4 TFZ=7.1 PAK=38 LLG=208 TFZ==7.3 RFZ=4.0 TFZ=8.5 PAK=38
    RFZ=3.2 TFZ=8.0 PAK=38 LLG=308 TFZ==8.5
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.3   62.2  219.9 FRAC -0.24  0.20 -0.33 BFAC  6.17
   SOLU 6DIM ENSE PDB_idealhelix EULER  275.9  114.0   38.4 FRAC -0.14  0.10 -0.05 BFAC  6.88 #TFZ==7.3
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.09 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  277.1   43.5  223.6 FRAC -0.04 -0.20 -0.34 BFAC  9.78 #TFZ==8.5
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

   Partial Solution #10 annotation (history):
   SOLU SET  PAK=0 RFZ=4.2 TFZ=5.0 PAK=0 LLG=142 RFZ=3.8 TFZ=7.4 PAK=39 LLG=211 TFZ==9.4 RFZ=4.0 TFZ=8.5 PAK=39 RFZ=3.0
    TFZ=6.7 PAK=39 LLG=307 TFZ==8.3
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE FIX EULER    0.0    0.0    0.0 FRAC -0.00 -0.00  0.00 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER   16.2  122.6   39.8 FRAC -0.16  0.07 -0.03 BFAC  6.98
   SOLU 6DIM ENSE PDB_idealhelix EULER   15.0   62.8  219.4 FRAC -0.25  0.20 -0.33 BFAC  5.62 #TFZ==9.4
   SOLU 6DIM ENSE PDB_idealhelix EULER  139.2  107.7   52.5 FRAC  0.10 -0.08 -0.11 BFAC  0.00
   SOLU 6DIM ENSE PDB_idealhelix EULER  215.3   51.2  243.3 FRAC -0.00 -0.24 -0.39 BFAC 10.13 #TFZ==8.3
   SOLU ENSEMBLE FIX VRMS DELTA +0.0010 #RMSD  0.10 #VRMS  0.11

CPU Time: 0 days 3 hrs 10 mins 25.60 secs (  11425.60 secs)
Finished: Thu Mar 12 04:02:50 2020

--------
ADVISORY
--------


The correlation corrected RmsD of ensemble "PDB_idealhelix" is 0.305 but the largest RmsD
assigned to models in the ensemble is only 0.100. Is this intentional?
eLLG indicates that placement of a single copy of ensemble "PDB_idealhelix" will be very
difficult


eLLG indicates that best placement of ensemble "PDB_idealhelix" will definitely be correct
in the context of already correctly placed components

The top solution from a FTF did not pack

The top solution from a TF rescoring did not pack

--------
WARNINGS
--------

------------------------------------------------------------------------------------------
Warning: Patterson Pathology in NMOL analysis, possible coiled coil or other regular
structure, view Patterson (TNCS PATT MAPS ON)
------------------------------------------------------------------------------------------

CPU Time: 0 days 3 hrs 10 mins 25.60 secs (  11425.60 secs)
Finished: Thu Mar 12 04:02:50 2020
"""
        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = PhaserParser(stdout=stdout_contents, fname=fname)
        parser.parse()

        self.assertEqual('318', parser.LLG)
        self.assertEqual('12.0', parser.eLLG)
        self.assertEqual('10.2', parser.TFZ)
        self.assertEqual('3.1', parser.RFZ)
        self.assertEqual('+0.0010', parser.VRMS)

        self.assertTupleEqual(('318', '10.2', '3.1', '12.0', '+0.0010'), parser.summary)


if __name__ == '__main__':
    unittest.main()
