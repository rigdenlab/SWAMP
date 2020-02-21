import os
import unittest
from swamp.wrappers.wrefmac import wRefmac


class MockRefmac(wRefmac):
    """A class to mock :py:obj:`~swmap.wrappres.wrefmac.wRefmac` for testing purposes"""

    def make_workdir(self):
        """Override :py:func:`~swmap.wrappres.wrefmac.wRefmac.make_workdir` so that it does not create \
        :py:attr:`~swmap.wrappres.wrefmac.wRefmac.workdir`"""
        pass

    def run(self):
        """Override :py:func:`~swmap.wrappres.wrefmac.wRefmac.run` so that it does not execute \
        :py:attr:`~swmap.wrappres.wrefmac.wRefmac.cmd`"""

        self.logcontents = """<B><FONT COLOR="#FF0000"><!--SUMMARY_BEGIN-->
<html> <!-- CCP4 HTML LOGFILE -->
<hr>
<!--SUMMARY_END--></FONT></B>
<B><FONT COLOR="#FF0000"><!--SUMMARY_BEGIN-->
<pre>
 
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.1.000: Refmac          version 5.8.0258 : 09/10/19##
 ###############################################################
 User: unknown  Run date: 17/ 2/2020 Run time: 02:55:53 


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.

 $$
 $SUMMARY :Reference4:  $$ Refmac: $$
 :TEXT:Reference2: $$

  Data line--- RIDG DIST SIGM  0.02
  Data line--- MAKE HYDR N
  Data line--- WEIGHT MATRIX 0.01
  Data line--- NCYC 100
  Data line--- END

 OPENED INPUT MTZ FILE 
 Logical Name: HKLIN   Filename: /data1/filo/results/SWAMP_benchmarking/3zux.mtz 

 LABIN FP=FP SIGFP=SIGFP FREE=FREE

    ****                     Input and Default parameters#                      ****


Input coordinate file.  Logical name - XYZIN actual file name  - /data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/phaser/phaser_out.pdb
Output coordinate file. Logical name - XYZOUT actual file name - /data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/refmac/refmac_out.pdb
Input reflection file.  Logical name - HKLIN actual file name  - /data1/filo/results/SWAMP_benchmarking/3zux.mtz
Output reflection file. Logical name - HKLOUT actual file name - /data1/filo/results/SWAMP_benchmarking/SWAMP_0/swamp_mr/search_461/run_1/refmac/refmac_out.mtz

Cell from mtz :    73.330    73.330   163.520    90.000    90.000    90.000
Space group from mtz: number -   91; name - P 41 2 2

  Refinement type                        : Restrained

    ****                           Makecif parameters                           ****

Dictionary files for restraints : /home/filo/opt/CCP4_7.1/ccp4-7.1/lib/data/monomers/mon*cif
Parameters for new entry and VDW: /home/filo/opt/CCP4_7.1/ccp4-7.1/lib/data/monomers/ener_lib.cif
    Apart from amino acids and DNA/RNA all monomers will be checked to see if atom names and connectivity is correct
    Hydrogens will not be used
    Links between monomers will be checked. Only those links present in the coordinate file will be used
    Standard sugar links will be analysed and used
    For new ligands "ideal restraint" values will be taken from the energetic libary ener_lib.cif
    Symmetry related links will be analysed and used
    Cis peptides will be found and used automatically


  Residual                               : Rice Maximum Likelihood for Fs

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
                        Standard  External       All
                Bonds:       208         0       208
               Angles:       288         0       288
              Chirals:        42         0        42
               Planes:        40         0        40
             Torsions:       162         0       162
            Intervals:         0         0         0
--------------------------------------------------------------------------------


  Number of harmonic restraints       =            0
  Number of atoms in special position =            0
 -----------------------------------------------------


    ****                           Info from mtz file                           ****

FreeR flag                            0
Number of "free" reflections       1180
Number of   all  reflections      23385
--------------------------------------------------------------------------------
 Number of reflections in file      23385
 Number of reflections read          23385

   Current auto weighting coefficient =    8.0001259    

 mode : HKRF


######  TLS Group Definitions ######

Resolution limits                    =     24.443     2.200
Number of used reflections           =      22205
Percentage observed                  =    99.5318
Percentage of free reflections       =     5.0460
Overall R factor                     =     0.5454
Free R factor                        =     0.5742
Average Fourier shell correlation    =     0.3949
AverageFree Fourier shell correlation=     0.3721
Overall weighted R factor            =     0.5453
Free weighted R factor               =     0.5780
Overall weighted R2 factor           =     0.6532
Free weighted R2 factor              =     0.7662
Average correlation coefficient      =     0.1512
Overall correlation coefficient      =     0.5718
Free correlation coefficient         =     0.4532
Cruickshanks DPI for coordinate error=     0.1212
DPI based on free R factor           =     0.1252
Overall figure of merit              =     0.3010
ML based su of positional parameters =     0.3464
ML based su of thermal parameters    =    12.7396
-----------------------------------------------------------------------------
  Time in seconds: CPU =       375.55
             Elapsed =         384.00

    ****           Things for loggraph, R factor and others vs cycle            ****


$TABLE: Rfactor analysis, stats vs cycle  :
$GRAPHS:<Rfactor> vs cycle :N:1,2,3:
:FOM vs cycle :N:1,4:
:-LL vs cycle :N:1,5:
:-LLfree vs cycle :N:1,6:
:Geometry vs cycle:N:1,7,8,9,10,11:
$$
    Ncyc    Rfact    Rfree     FOM      -LL     -LLfree  rmsBOND  zBOND rmsANGL  zANGL rmsCHIRAL $$
$$
       0   0.5655   0.5811   0.292      108413.    5720.6   0.0127  0.983   1.385  0.906   0.037
       1   0.5587   0.5789   0.286      108253.    5722.3   0.0069  0.592   1.122  0.733   0.060
       2   0.5568   0.5793   0.282      108126.    5718.4   0.0065  0.570   1.161  0.758   0.067
 $$
 $TEXT:Result: $$ Final results $$
                      Initial    Final
           R factor    0.5655   0.5454
             R free    0.5811   0.5742
     Rms BondLength    0.0127   0.0067
      Rms BondAngle    1.3854   1.5617
     Rms ChirVolume    0.0365   0.1007
 $$
 Harvest: NO PNAME_KEYWRD given - no deposit file created
<B><FONT COLOR="#FF0000"><!--SUMMARY_BEGIN-->
 Refmac:  End of Refmac_5.8.0258   
Times: User:     390.4s System:    1.1s Elapsed:     6:32  
</pre>
</html>
<!--SUMMARY_END--></FONT></B>
"""
        self.get_scores()


class MyTestCase(unittest.TestCase):

    def test_1(self):
        refmac = MockRefmac(mtzin='/empty/path/fname.mtz', pdbin='/empty/path/fname.pdb', ridg_dist_sigm="0.52",
                            workdir='/empty/path/workdir', make_hydr='Y', weight_matrix="0.21", ncyc="140")

        self.assertListEqual(refmac.cmd, ["refmac5", 'hklin', '/empty/path/fname.mtz', 'hklout',
                                          '/empty/path/workdir/refmac_out.mtz', 'xyzin', '/empty/path/fname.pdb',
                                          'xyzout', '/empty/path/workdir/refmac_out.pdb'])

        self.assertEqual("RIDG DIST SIGM  0.52" + os.linesep + "MAKE HYDR Y" + os.linesep + "WEIGHT MATRIX 0.21" \
                         + os.linesep + "NCYC 140" + os.linesep + "END", refmac.keywords)

        refmac.run()

        self.assertTupleEqual(('0.0365', '0.1007'), refmac.chirvol_delta)
        self.assertTupleEqual(('1.3854', '1.5617'), refmac.bondangle_delta)
        self.assertTupleEqual(('0.0127', '0.0067'), refmac.bondlenght_delta)
        self.assertTupleEqual(('0.5655', '0.5454'), refmac.rfactor_delta)
        self.assertTupleEqual(('0.5811', '0.5742'), refmac.rfree_delta)
        self.assertEqual('0.5742', refmac.rfree)
        self.assertEqual('0.5454', refmac.rfactor)
        self.assertEqual("Refmac results: Rfactor - 0.5454   Rfree - 0.5742   Local CC - NA   Overall CC - NA",
                         refmac.summary_results)


if __name__ == '__main__':
    unittest.main()
