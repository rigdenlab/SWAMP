import os
import unittest
import numpy as np
import pandas as pd
from swamp.utils import create_tempfile
from swamp.wrappers.gesamt import Gesamt


class MockGesamt(Gesamt):
    """A class to mock :py:obj:`~swmap.wrappers.gesamt.Gesamt` for testing purposes"""

    def run(self):
        """Override :py:func:`~swmap.wrappers.gesamt.Gesamt.run` for testing purposes"""

        if self.mode == "search-archive":
            file_contents = """#  Hit   PDB  Chain  Q-score  r.m.s.d     Seq.  Nalign  nRes    File
#  No.   code   Id                         Id.                  name
     1          A     1.0000   0.0000   1.0000     47     47   6qd5_5A_17A.pdb
     2          A     0.4736   0.0000   1.0000     32     46   6qd5_5A_15A.pdb
     3          A     0.4736   0.0000   1.0000     32     46   6qd5_5A_9A.pdb
"""
            self.hits_out = create_tempfile(content=file_contents)

        elif self.mode == "alignment":
            if len(self.pdbin) == 2:

                self.logcontents = """
 GESAMT: General Efficient Structural Alignment of Macromolecular Targets
 ------------------------------------------------------------------------
 Version 1.15 of 25-Jan-2017, built with MMDB v.2.0.17
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.073: Gesamt           version 7.0.073 :         ##
 ###############################################################
 User: filo  Run date: 22/ 2/2020 Run time: 13:12:47 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

$TEXT:Reference: $$Please cite$$
E. Krissinel (2012). Enhanced fold recognition using efficient
short fragment clustering. J. Mol. Biochem. 1(2) 76-85.
$$
<!--SUMMARY_BEGIN-->

 ===============================================================================

 ... reading FIXED structure : file '/mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb', selection '*'
        308 atoms selected
 ... reading MOVING structure: file '/mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb', selection '*'
        355 atoms selected
<!--SUMMARY_END-->

 ===============================================================================

 CPU stage 1 (clustering):   0.06758 secs
 CPU stage 2 (refinement):   0.02398 secs

 ===== Structures

     Ref.  |  Nres  | File (selection)                            
   ========+========+=============================================
     FIXED |   308  | /mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb (*) 
    MOVING |   355  | /mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb (*) 

 have been aligned and superposed.


 ===============================================================================

 SUPERPOSITION
 ~~~~~~~~~~~~~

 Q-score          : 0.029     
 RMSD             : 3.357     
 Aligned residues : 84
 Sequence Id:     : 0.048     

 Transformation matrix for FIXED structure is identity.

 Transformation matrix for MOVING structure:

          Rx           Ry           Rz             T
      0.95978      0.09955      0.26253      -15.08154
     -0.23123     -0.25016      0.94019       84.32709
      0.15927     -0.96307     -0.21708       83.61516

 Direction cosines of the rotation axis: -0.98383 0.05338 -0.17098
 Rotation angle                        : 104.69832 

 in fractional coordinates of FIXED structure:

          Rx           Ry           Rz             T
      0.95978      0.09955      0.58542       -0.20567
     -0.23123     -0.25016      2.09655        1.14997
      0.07142     -0.43189     -0.21708        0.51135

 in fractional coordinates of MOVING structure:

          Rx           Ry           Rz             T
      0.95978      0.16017      0.27547       -0.16597
     -0.14371     -0.25016      0.61314        0.57675
      0.15178     -1.47678     -0.21708        0.87693

 -------------------------------------------------------------------------------

 CENTROIDS
 ~~~~~~~~~            Orthogonal                       Fractional
               X          Y          Z            XF       YF       ZF
  FIXED     18.19964   69.74116   81.11004      0.24819  0.95106  0.49603
 MOVING     18.75503   20.49314   -5.04154      0.20639  0.14016 -0.05287

 Distance between centroids                   : 99.23592  
 Direction cosines of vector between centroids: 0.00560 -0.49627 -0.86815
 Angle between rotation axis and vector between centroids: 83.31328 

 -------------------------------------------------------------------------------

 CCP4 format rotation-translation operator
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 Polar angles (omega,phi,kappa)   :     99.84487    176.89438    104.69832
 Euler angles (alpha,beta,gamma)  :     74.39866    102.53743    -99.39011
 Orthogonal translation (Angstrom):    -15.08154     84.32709     83.61516

 ===============================================================================

 RESIDUE ALIGNMENT
 ~~~~~~~~~~~~~~~~~
$$

.-------------.------------.-------------.
|   FIXED     |  Dist.(A)  |   MOVING    |
|-------------+------------+-------------|
|             |            | - A:PHE  -2 |
|             |            | + A:GLN  -1 |
|             |            | . A:SER   0 |
|H- A:ALA 307 |            |             |
|H+ A:LYS 308 |            |             |
| - A:ALA 309 |            |             |
`-------------'------------'-------------'

 Notations:
 S/H   residue belongs to a strand/helix
 +/-/. hydrophylic/hydrophobic/neutral residue
 **    identical residues matched: similarity 5
 ++    similarity 4
 ==    similarity 3
 --    similarity 2
 ::    similarity 1
 ..    dissimilar residues: similarity 0
 Gesamt:  Normal termination
"""
            else:
                self.logcontents = """
 GESAMT: General Efficient Structural Alignment of Macromolecular Targets
 ------------------------------------------------------------------------
 Version 1.15 of 25-Jan-2017, built with MMDB v.2.0.17
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.073: Gesamt           version 7.0.073 :         ##
 ###############################################################
 User: filo  Run date: 22/ 2/2020 Run time: 13:14:47 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

$TEXT:Reference: $$Please cite$$
E. Krissinel (2012). Enhanced fold recognition using efficient
short fragment clustering. J. Mol. Biochem. 1(2) 76-85.
$$
<!--SUMMARY_BEGIN-->

 ===========================================================
 ... reading file '/mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb', selection '*':
        308 atoms selected
 ... reading file '/mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb', selection '*':
        355 atoms selected
 ... reading file '/mnt/sdb1/SWAMP_benchmark/4njn/4njn.pdb', selection '*':
        182 atoms selected

 Parameter Q-score:                3.000 angstroms
 Weighted superposition is not used
 Number of threads used:           1


 CPU stage 1 (cross-alignments):   0.20074 secs
 CPU stage 2 (refinement):         0.00997 secs


 ===== Structures

     Ref.  |  Nres  | File (selection)                            
   ========+========+=============================================
      S001 |   308  | /mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb (*) 
      S002 |   355  | /mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb (*) 
      S003 |   182  | /mnt/sdb1/SWAMP_benchmark/4njn/4njn.pdb (*) 

 have been aligned and superposed.


 ===== Superposition matrices (orthogonal):

   ____________________________________________________________________
   (o) For structure S001 [/mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb(*)]:

        Rx         Ry         Rz           T
      1.000     -0.000     -0.000       -0.000
      0.000      1.000      0.000        0.000
      0.000      0.000      1.000        0.000

   ____________________________________________________________________
   (o) For structure S002 [/mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb(*)]:

        Rx         Ry         Rz           T
      0.817     -0.080     -0.571      -24.547
      0.575      0.189      0.796       64.490
      0.044     -0.979      0.200       86.418

   ____________________________________________________________________
   (o) For structure S003 [/mnt/sdb1/SWAMP_benchmark/4njn/4njn.pdb(*)]:

        Rx         Ry         Rz           T
      0.810     -0.569     -0.144       31.572
      0.510      0.803     -0.309       76.523
      0.292      0.177      0.940       34.279

 ===== Superposition matrices (fractional):

 ===== Scores achieved:

   quality Q:  0.0099  (normalised to [0...1])
     r.m.s.d:  2.1311  (A)
      Nalign:  44      (residues)

   ______________________________________________________
   (o) pairwise Q-scores (consensus Q-score on diagonal):

             S001    S002    S003 
         .------------------------
     S001|  0.045   0.002   0.006  
     S002|  0.002   0.036   0.005  
     S003|  0.006   0.005   0.101  

   _______________________________________________________
   (o) pairwise r.m.s.d. (consensus r.m.s.d. on diagonal):

             S001    S002    S003 
         .------------------------
     S001|  2.212   4.225   3.234  
     S002|  4.225   2.366   3.545  
     S003|  3.234   3.545   1.770  

   _____________________
   (o) pairwise seq. Id:

             S001    S002    S003 
         .------------------------
     S001|  1.000   0.091   0.091  
     S002|  0.091   1.000   0.000  
     S003|  0.091   0.000   1.000  


 ===== Residue alignment:

  Disp. | |    S001    | |    S002    | |    S003    
 -------+-+------------+-+------------+-+------------
        | |            | |  A:PHE  -2 | |            
        | |            | |  A:GLN  -1 | |            
        | |            | |  A:SER   0 | |            
        | |            | |  A:MET   1 | |            
        | |  A:ALA 309 | |            | |            
 -------'-'------------'-'------------'-'------------
 Gesamt:  Normal termination
"""
        self.get_scores()


class GesamtWrapperTestCase(unittest.TestCase):

    def test_1(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='alignment', pdbout='/empty/path/pdbout',
                            pdbin=('/empty/path/pdb_a', '/empty/path/pdb_b'))

        self.assertListEqual(gesamt.cmd, [os.path.join(os.environ['CCP4'], 'bin', 'gesamt'), '/empty/path/pdb_a',
                                          '/empty/path/pdb_b', '-o', '/empty/path/pdbout'])
        gesamt.run()
        self.assertTupleEqual((gesamt.qscore, gesamt.rmsd, gesamt.seq_id, gesamt.n_align), (0.029, 3.357, 0.048, 84))

    def test_2(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='alignment', pdbout='/empty/path/pdbout',
                            pdbin=('/empty/path/pdb_a', '/empty/path/pdb_b', '/empty/path/pdb_c'))

        self.assertListEqual(gesamt.cmd, [os.path.join(os.environ['CCP4'], 'bin', 'gesamt'), '/empty/path/pdb_a',
                                          '/empty/path/pdb_b', '/empty/path/pdb_c', '-o', '/empty/path/pdbout'])
        gesamt.run()
        self.assertTupleEqual((gesamt.qscore, gesamt.rmsd, gesamt.seq_id, gesamt.n_align), (0.0099, 2.1311, np.nan, 44))

    def test_3(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='alignment', pdbout='/empty/path/pdbout',
                            pdbin=('/empty/path/pdb_a',))

        self.assertIsNone(gesamt.cmd)
        self.assertTrue(gesamt.error)

    def test_4(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='search-archive', pdbin='/empty/path/pdbin', min2=0.7,
                            gesamt_archive='/empty/path/gesamt-archive', hits_out='/empty/path/out.hits', nthreads=15,
                            min1=0.4)

        self.assertListEqual(gesamt.cmd, [os.path.join(os.environ['CCP4'], 'bin', 'gesamt'), '/empty/path/pdbin',
                                          '-archive', '/empty/path/gesamt-archive', '-o', '/empty/path/out.hits',
                                          '-nthreads=15', '-min1=0.4', '-min2=0.7'])
        self.assertListEqual(gesamt._filthy_files, ['/empty/path/out.hits'])
        gesamt.run()
        df = pd.DataFrame([['1.0000', '0.0000', '1.0000', '47', '47', '6qd5_5A_17A.pdb'],
                           ['0.4736', '0.0000', '1.0000', '32', '46', '6qd5_5A_15A.pdb'],
                           ['0.4736', '0.0000', '1.0000', '32', '46', '6qd5_5A_9A.pdb']])
        df.columns = ["qscore", "rmsd", "seq_id", "n_align", "n_res", "fname"]

        try:
            pd.testing.assert_frame_equal(df, gesamt.summary_results)
        except AssertionError as e:
            raise self.failureException('Non identical pd.DataFrame') from e

        self.addCleanup(os.remove, gesamt.hits_out)

    def test_5(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='search-archive', pdbin=None, min2=0.7, nthreads=15,
                            gesamt_archive='/empty/path/gesamt-archive', hits_out='/empty/path/out.hits', min1=0.4)
        self.assertTrue(gesamt.error)

    def test_6(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='unknown', pdbin=None, min2=0.7, nthreads=15,
                            gesamt_archive='/empty/path/gesamt-archive', hits_out='/empty/path/out.hits', min1=0.4)
        self.assertTrue(gesamt.error)

    def test_7(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='make-archive', pdb_archive='/empty/path/pdb-archive',
                            gesamt_archive='/empty/path/gesamt-archive')

        self.assertListEqual(gesamt.cmd, [os.path.join(os.environ['CCP4'], 'bin', 'gesamt'), '--make-archive',
                                          '/empty/path/gesamt-archive', '-pdb', '/empty/path/pdb-archive'])

    def test_8(self):
        gesamt = MockGesamt(workdir='/empty/path/workdir', mode='make-archive', pdb_archive=None,
                            gesamt_archive='/empty/path/gesamt-archive')

        self.assertIsNone(gesamt.cmd)
        self.assertTrue(gesamt.error)
