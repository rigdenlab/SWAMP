import unittest
from swamp.parsers.gesamtparser import GesamtParser, GesamtErrorCodes


class GesamtParserTestCase(unittest.TestCase):

    def test_1(self):
        STDOUT = """
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
        result = GesamtParser.get_pairwise_qscores(STDOUT)
        self.assertDictEqual(result, {'/mnt/sdb1/SWAMP_benchmark/3zux/3zux.pdb': 0.045,
                                      '/mnt/sdb1/SWAMP_benchmark/4njn/4njn.pdb': 0.101,
                                      '/mnt/sdb1/SWAMP_benchmark/5mlz/5mlz.pdb': 0.036})

    def test_2(self):
        STDOUT = """ GESAMT: General Efficient Structural Alignment of Macromolecular Targets
 ------------------------------------------------------------------------
 Version 1.16 of 14-Jan-2020, built with MMDB v.2.0.20
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.1.000: Gesamt           version 7.1.000 :         ##
 ###############################################################
 User: filo  Run date: 23/ 2/2020 Run time: 21:28:28 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

$TEXT:Reference: $$Please cite$$
E. Krissinel (2012). Enhanced fold recognition using efficient
short fragment clustering. J. Mol. Biochem. 1(2) 76-85.
$$
<!--SUMMARY_BEGIN-->

 ===============================================================================

 ... reading FIXED structure : file '/home/filo/opt/CCP4/ccp4-7.0/lib/py2/swamp/library/db/pdb/5fvn_27A_29A.pdb', selection '*'
          0 atoms selected with warning (rc=15)
      ***** ERROR #15 READ:

 File can not be opened.


            crystal data not found
 ... reading MOVING structure: file '/home/filo/opt/CCP4/ccp4-7.0/lib/py2/swamp/library/db/pdb/5fvn_23A_25A.pdb', selection '*'
          0 atoms selected with warning (rc=15)
      ***** ERROR #15 READ:

 File can not be opened.


            crystal data not found
<!--SUMMARY_END-->

 ===============================================================================

 CPU stage 1 (clustering):   0.00001 secs
 CPU stage 2 (refinement):   0.00000 secs

 ===============================================================================

 Query      /home/filo/opt/CCP4/ccp4-7.0/lib/py2/swamp/library/db/pdb/5fvn_27A_29A.pdb(*)
 and Target /home/filo/opt/CCP4/ccp4-7.0/lib/py2/swamp/library/db/pdb/5fvn_23A_25A.pdb(*)

 are DISSIMILAR and cannot be reasonably aligned.

 ===============================================================================
 Gesamt:  Normal termination
"""
        parser = GesamtParser(mode='alignment', stdout=STDOUT)
        parser.parse()
        self.assertTrue(parser.error)
        self.assertEqual(parser.error, GesamtErrorCodes.DISSIMILAR)

    def test_3(self):
        STDOUT = """
 GESAMT: General Efficient Structural Alignment of Macromolecular Targets
 ------------------------------------------------------------------------
 Version 1.16 of 14-Jan-2020, built with MMDB v.2.0.20
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.1.000: Gesamt           version 7.1.000 :         ##
 ###############################################################
 User: filo  Run date: 23/ 2/2020 Run time: 21:43:27 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

$TEXT:Reference: $$Please cite$$
E. Krissinel (2012). Enhanced fold recognition using efficient
short fragment clustering. J. Mol. Biochem. 1(2) 76-85.
$$
<!--SUMMARY_BEGIN-->

 ===========================================================
 ... reading file '/home/filo/opt/CCP4/ccp4-7.0/lib/py2/swamp/library/db/pdb/5fvn_29A_31A.pdb', selection '*':
          0 atoms selected with warning (rc=15)
      ***** ERROR #15 READ:

 File can not be opened.



 Parameter Q-score:                3.000 angstroms
 Weighted superposition is not used
 Number of threads used:           1



 STOP DUE TO READ ERRORS
 --- check input file format

 ===========================================================
 Gesamt:  Normal termination

"""
        parser = GesamtParser(mode='alignment', stdout=STDOUT)
        parser.parse()
        self.assertTrue(parser.error)
        self.assertEqual(parser.error, GesamtErrorCodes.READ_ERRORS)

    def test_4(self):
        STDOUT = """
 GESAMT: General Efficient Structural Alignment of Macromolecular Targets
 ------------------------------------------------------------------------
 Version 1.15 of 25-Jan-2017, built with MMDB v.2.0.17
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.073: Gesamt           version 7.0.073 :         ##
 ###############################################################
 User: filo  Run date: 23/ 2/2020 Run time: 21:49:18 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

$TEXT:Reference: $$Please cite$$
E. Krissinel (2012). Enhanced fold recognition using efficient
short fragment clustering. J. Mol. Biochem. 1(2) 76-85.
$$
<!--SUMMARY_BEGIN-->

 ===========================================================
 ... reading file '/home/filo/opt/CCP4/ccp4-7.0/share/swamp/pdb/5fvn_29A_31A.pdb', selection '*':
         17 atoms selected
 ... reading file '/home/filo/opt/CCP4/ccp4-7.0/share/swamp/pdb/5fvn_19A_21A.pdb', selection '*':
         17 atoms selected
 ... reading file '/home/filo/opt/CCP4/ccp4-7.0/share/swamp/pdb/5fvn_17A_19A.pdb', selection '*':
         17 atoms selected
 ... reading file '/home/filo/opt/CCP4/ccp4-7.0/share/swamp/pdb/5fvn_15A_17A.pdb', selection '*':
         18 atoms selected

 Parameter Q-score:                3.000 angstroms
 Weighted superposition is not used
 Number of threads used:           1


 CPU stage 1 (cross-alignments):   0.00005 secs
 CPU stage 2 (refinement):         0.00000 secs



 ALIGNMENT ERROR 2
 ===========================================================
 Gesamt:  Normal termination
"""
        parser = GesamtParser(mode='alignment', stdout=STDOUT)
        parser.parse()
        self.assertTrue(parser.error)
        self.assertEqual(parser.error, GesamtErrorCodes.ERROR_2)


if __name__ == '__main__':
    unittest.main()
