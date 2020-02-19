import os
import unittest
from swamp.utils import create_tempfile
from swamp.parsers.shelxeparser import ShelxeParser


class MyTestCase(unittest.TestCase):

    def test_1(self):
        stdout_contents = """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     +  SHELXE  -  PHASING AND DENSITY MODIFICATION  -  Version 2019/1  +
     +  Copyright (c)  George M. Sheldrick and Isabel Uson 2001-19      +
     +  Started at 03:02:43 on 17 Feb 2020                              +
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     Please cite: I. Uson & G.M. Sheldrick (2018), "An introduction to
     experimental phasing of macromolecules illustrated by SHELX;
     new autotracing features" Acta Cryst. D74, 106-116
     (Open Access) if SHELXE proves useful.

     Command line parameters:
      shelxe-input.pda -s0.6 -a20 -m10 -t4 -l2.0 -n -q -f -o

     Cell and symmetry from shelxe-input.pda
     Phases calculated using atoms from shelxe-input.pda
     Native data from shelxe-input.hkl
     Pruned PDB fragment output to shelxe-input.pdo
     Listing output to shelxe-input.lst
     Phases output to shelxe-input.phs
     Poly-Ala trace output to shelxe-input.pdb

     Summary of parameters to be employed:

     -a    20  global autotracing cycles
     -b   5.0  extra B for revised heavy atom sites
     -B unset  just build single beta-strands
     -c 0.400  fraction of pixels in crossover region
     -d 0.000  high resolution limit to be applied to input data
     -D unset  do not fuse disulfides for NCS
     -e unset  fill in missing data up to maximum resolution + 0.2 Ang.
     -f        read F instead of intensity from native .hkl file
     -F 0.800  fractional weight for phases from previous global cycle
     -g 1.100  solvent gamma flipping factor
     -G 0.700  FOM threshold for initial tripeptides and chain building
     -i unset  no structure inversion
     -k   4.5  minimum height/sigma for revised heavy atom sites
     -l     2  space for  2000000 reflections
     -L     6  minimum number of residues/chain (if more than 3 chains)
     -m    10  cycles of density modification
     -o        omit residues from fragment to optimize CC
     -p unset  no phosphate search
     -q     7   alpha-helix search
     -r  3.00  map resolution (multiplies maximum indices)
     -s 0.600  solvent fraction
     -S  2.42  radius of sphere of influence
     -t  4.00  time for initial searches (-t6 or more if difficult)
     -u   500  MB allocatable memory for fragment optimization
     -U  0.00  abort if less than this % of fragment CA retained within 0.7A
     -v 0.000  density sharpening factor
     -w 0.200  weight for experimental phases after cycle 1
     -x unset  no phase and trace diagnostics 
     -y  1.80  highest resolution in Ang. for starting phases from model
     -z unset  do not optimize heavy atoms

     Space group: P 41 2 2   Allowed origin shift code:  9

        23385 Reflections read from file shelxe-input.hkl

        23385 Unique data, highest resolution =  2.200 Angstroms

     Anisotropic scaling: intensities multiplied by
      0.000624h^2 +0.000624k^2 -0.000100l^2 +0.000000kl +0.000000hl +0.000000hk

        72 Reflections with d > 2.400 and      0 in range  2.400 > d > 2.200 added
     Fourier grid =   128 x  128 x   29     0.000 <= z <= 0.125

        92 Point spherical net set up with radius 2.42A
        24 Extra Fourier layers will be generated


       210 Atoms read from PDB format file shelxe-input.pda

         0 NCS groups defined in .pda input file
      17 Residue blocks eliminated out of  42, increasing CC from  8.82 to 11.40%
      25 residues written to shelxe-input.pdo

     Overall CC between native Eobs and Ecalc (from fragment) = 13.93%
     Starting from fragment phases truncated to  1.800A
     <wt> = 0.295 for fragment phases
     <wt> = 0.295, Contrast = 0.196, Connect. = 0.617 for dens.mod. cycle 1
     <wt> = 0.296, Contrast = 0.244, Connect. = 0.628 for dens.mod. cycle 2
     <wt> = 0.297, Contrast = 0.265, Connect. = 0.642 for dens.mod. cycle 3
     <wt> = 0.299, Contrast = 0.279, Connect. = 0.652 for dens.mod. cycle 4
     <wt> = 0.300, Contrast = 0.289, Connect. = 0.659 for dens.mod. cycle 5
     <wt> = 0.300, Contrast = 0.298, Connect. = 0.665 for dens.mod. cycle 6
     <wt> = 0.300, Contrast = 0.305, Connect. = 0.671 for dens.mod. cycle 7
     <wt> = 0.300, Contrast = 0.311, Connect. = 0.676 for dens.mod. cycle 8
     <wt> = 0.300, Contrast = 0.316, Connect. = 0.679 for dens.mod. cycle 9
     <wt> = 0.300, Contrast = 0.320, Connect. = 0.682 for dens.mod. cycle 10

     NOGO map generated for regions about rotation axes (if any)

      1345 peaks > 0.5 sigma used to seed fragment search
     Space for about  350 unique residues taking solvent into account

       139 potential (overlapping) alpha-helix fragments found

       #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  

      56.0% starting CA preserved within 1.0A,  52.0% / 0.7A and  44.0% / 0.5A

     CC for partial structure against native data =  14.46 %

     ------------------------------------------------------------------------------

     Global autotracing cycle   2

     <wt> = 0.468 for fragment phases
     <wt> = 0.300, Contrast = 0.246, Connect. = 0.649 for dens.mod. cycle 1
     <wt> = 0.300, Contrast = 0.274, Connect. = 0.651 for dens.mod. cycle 2
     <wt> = 0.300, Contrast = 0.310, Connect. = 0.679 for dens.mod. cycle 3
     <wt> = 0.300, Contrast = 0.324, Connect. = 0.689 for dens.mod. cycle 4
     <wt> = 0.300, Contrast = 0.333, Connect. = 0.694 for dens.mod. cycle 5
     <wt> = 0.300, Contrast = 0.338, Connect. = 0.698 for dens.mod. cycle 6
     <wt> = 0.300, Contrast = 0.341, Connect. = 0.700 for dens.mod. cycle 7
     <wt> = 0.300, Contrast = 0.345, Connect. = 0.702 for dens.mod. cycle 8
     <wt> = 0.300, Contrast = 0.347, Connect. = 0.702 for dens.mod. cycle 9
     <wt> = 0.300, Contrast = 0.348, Connect. = 0.703 for dens.mod. cycle 10

     NOGO map generated for regions about rotation axes (if any)

      1331 peaks > 0.5 sigma used to seed fragment search
     Space for about  350 unique residues taking solvent into account

       206 potential (overlapping) alpha-helix fragments found

       #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  
     A:   13  24.554   2.434   0.819   0.545   0.833   0.500   CB 1.839
     B:   12  16.211   2.379   0.435   0.486   0.818   0.636   CB 1.775
     C:    9  12.906   1.744   0.682   0.769   0.875   0.500   CB 1.400
     D:   10   5.099   1.504   0.055   0.682   0.778   0.667   CB 1.165
     E:    8   4.370   1.356   0.103   0.746   0.857   0.571   CB 1.034
     F:    8   6.876   1.402   0.601   0.630   0.857   0.571   CB 1.196
     G:    7   8.645   1.489   0.359   0.869   1.000   0.833   CB 1.274
     H:   10   6.926   1.398   0.612   0.656   0.778   0.667   CB 1.139
     I:    9   6.779   1.317   0.545   0.806   0.875   0.875   CB 1.031
     J:    8   3.589   1.222   0.425   0.770   0.857   0.714   CB 0.716
     K:    6   4.672   1.162   0.427   0.900   1.000   1.000   O 0.900
     L:    7   7.308   1.237   0.956   0.910   1.000   0.833   CB 0.875
     M:    8   2.651   1.108   0.184   0.652   0.714   0.714   N 0.929
     N:    6   3.770   1.071   0.268   0.895   1.000   1.000   N 0.891
     O:    8   4.202   1.158   0.376   0.759   0.857   0.714   CB 0.924
     P:    6   3.146   0.949   0.332   0.895   1.000   1.000   C 0.798
     Q:    6   2.465   0.990   0.161   0.887   1.000   1.000   CB 0.691

       104 potential tripeptides employed

       #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  
     R:   10   9.277   1.900   1.149   0.233   0.667   0.444   CB 1.550
     S:   22   6.717   1.791   0.549   0.149   0.524   0.238   CB 1.554
     T:    8   3.999   1.293   0.145   0.453   1.000   0.429   CB 1.071
           6   1.967   1.237   0.178   0.408   0.800   0.600   O 0.830 ?
     U:    8   3.333   1.624   0.353   0.284   0.857   0.429   CB 0.853
     V:   10   3.215   1.755   0.403   0.205   0.444   0.444   CB 1.390
     W:   14  12.759   2.349   0.371   0.378   0.769   0.385   CB 1.627
        16 residues pruned to eliminate duplicates
     X:    6   2.283   1.494  -0.022   0.347   0.800   0.400   CB 1.032
           6   1.950   1.143   0.400   0.379   0.800   0.600   CB 0.775 ?

     Using tripeptides from previous cycle as seeds

       141 residues left after pruning, divided into chains as follows:
     A:  12   B:   9   C:   8   D:   8   E:   7   F:  10   G:   9   H:   8   I:   7
     J:   8   K:   6   L:   6   M:   6   N:   6   O:  12   P:   8   Q:  11

      60.0% starting CA preserved within 1.0A,  56.0% / 0.7A and  32.0% / 0.5A

     CC for partial structure against native data =  17.52 %

      ------------------------------------------------------------------------------

     Global autotracing cycle  3

     Phases from autotracing cycle  2 used as input for final density modification

     <wt> = 0.538 for fragment phases
     <wt> = 0.300, Contrast = 0.239, Connect. = 0.665 for dens.mod. cycle 1
     <wt> = 0.300, Contrast = 0.268, Connect. = 0.668 for dens.mod. cycle 2
     <wt> = 0.300, Contrast = 0.314, Connect. = 0.708 for dens.mod. cycle 3
     <wt> = 0.300, Contrast = 0.331, Connect. = 0.720 for dens.mod. cycle 4
     <wt> = 0.300, Contrast = 0.341, Connect. = 0.722 for dens.mod. cycle 5
     <wt> = 0.300, Contrast = 0.344, Connect. = 0.723 for dens.mod. cycle 6
     <wt> = 0.300, Contrast = 0.349, Connect. = 0.724 for dens.mod. cycle 7
     <wt> = 0.300, Contrast = 0.350, Connect. = 0.724 for dens.mod. cycle 8
     <wt> = 0.300, Contrast = 0.353, Connect. = 0.723 for dens.mod. cycle 9
     <wt> = 0.300, Contrast = 0.355, Connect. = 0.723 for dens.mod. cycle 10

     Estimated mean FOM and mapCC as a function of resolution
     d    inf - 4.85 - 3.82 - 3.33 - 3.02 - 2.79 - 2.63 - 2.49 - 2.38 - 2.29 - 2.21
     <FOM>   0.526  0.502  0.439  0.406  0.358  0.286  0.310  0.339  0.315  0.265
     <mapCC> 0.942  0.863  0.765  0.753  0.674  0.547  0.593  0.637  0.480  0.389
     N        2347   2365   2321   2332   2448   2226   2452   2336   2259   2299

     Estimated mean FOM = 0.375   Pseudo-free CC = 43.82 %

     Best trace (cycle 3 with CC 20.64%) was saved as shelxe-input.pdb

     ==============================================================================

           CPU times required in seconds
           -----------------------------
             2.5 - Setup, data input and phasing
            25.3 - FFTs and peak-searches
            64.1 - Sphere of influence
             2.0 - Rest of density modification
          2515.9 - Alpha-helix search
          1218.5 - Tripeptide search
          1418.3 - Chain tracing
             0.0 - NCS analysis
           103.2 - B-value refinement for trace
             2.0 - Rest of tracing
             0.0 - Comparison with known structure


     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     +  SHELXE finished at 04:31:58      Total time:      5351.91 secs  +
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     """
        file_contents = """TITLE  shelxe-input.pdb   Cycle  19   CC = 20.64%   183  residues in  14  chains
    CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2
    SCALE1      0.013637  0.000000  0.000000        0.00000
    SCALE2      0.000000  0.013637  0.000000        0.00000
    SCALE3      0.000000  0.000000  0.006115        0.00000
    ATOM      1  CA  ALA A   1      24.175  11.677 -28.178  1.00 20.00           C
    ATOM      2  C   ALA A   1      22.865  12.325 -27.744  1.00 20.00           C
    ATOM      3  O   ALA A   1      21.841  12.145 -28.400  1.00 20.00           O
    ATOM      4  N   ALA A   2      22.893  13.079 -26.641  1.00 20.00           N
    ATOM      5  CA  ALA A   2      21.697  13.745 -26.137  1.00 20.00           C
    ATOM      6  CB  ALA A   2      21.427  15.034 -26.920  1.00 20.00           C
    ATOM      7  C   ALA A   2      20.456  12.858 -26.179  1.00 20.00           C
    ATOM      8  O   ALA A   2      19.348  13.357 -26.358  1.00 20.00           O
    ATOM      9  N   ALA A   3      20.639  11.545 -26.010  1.00 20.00           N
    ATOM     10  CA  ALA A   3      19.522  10.607 -26.031  1.00 20.00           C
    ATOM    909  N   ALA N   6      12.748  28.352  30.277  1.00 20.00           N
    ATOM    910  CA  ALA N   6      13.022  27.069  30.915  1.00 20.00           C
    ATOM    911  CB  ALA N   6      14.129  26.296  30.191  1.00 20.00           C
    ATOM    912  C   ALA N   6      13.407  27.321  32.369  1.00 20.00           C
    ATOM    913  O   ALA N   6      12.957  26.603  33.259  1.00 20.00           O
    ATOM    914  N   ALA N   7      14.237  28.339  32.613  1.00 20.00           N
    ATOM    915  CA  ALA N   7      14.668  28.669  33.968  1.00 20.00           C
    END 
    """

        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = ShelxeParser(fname=fname, stdout=stdout_contents)
        parser.parse()

        self.assertEqual("20.64", parser.cc)
        self.assertEqual("13.0", parser.acl)
        self.assertEqual("NO", parser.solution)
        self.assertEqual(3.0599999999999987, parser.average_cc_delta)
        self.assertEqual('13.93', parser.cc_eobs_ecalc)
        self.assertTupleEqual(('20.64', '13.0', '13.93', 3.0599999999999987, 'NO'), parser.summary)

    def test_2(self):
        stdout_contents = """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +  SHELXE  -  PHASING AND DENSITY MODIFICATION  -  Version 2019/1  +
 +  Copyright (c)  George M. Sheldrick and Isabel Uson 2001-19      +
 +  Started at 03:02:43 on 17 Feb 2020                              +
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 Please cite: I. Uson & G.M. Sheldrick (2018), "An introduction to
 experimental phasing of macromolecules illustrated by SHELX;
 new autotracing features" Acta Cryst. D74, 106-116
 (Open Access) if SHELXE proves useful.

 Command line parameters:
  shelxe-input.pda -s0.6 -a20 -m10 -t4 -l2.0 -n -q -f -o

 Cell and symmetry from shelxe-input.pda
 Phases calculated using atoms from shelxe-input.pda
 Native data from shelxe-input.hkl
 Pruned PDB fragment output to shelxe-input.pdo
 Listing output to shelxe-input.lst
 Phases output to shelxe-input.phs
 Poly-Ala trace output to shelxe-input.pdb

 Summary of parameters to be employed:

 -a    20  global autotracing cycles
 -b   5.0  extra B for revised heavy atom sites
 -B unset  just build single beta-strands
 -c 0.400  fraction of pixels in crossover region
 -d 0.000  high resolution limit to be applied to input data
 -D unset  do not fuse disulfides for NCS
 -e unset  fill in missing data up to maximum resolution + 0.2 Ang.
 -f        read F instead of intensity from native .hkl file
 -F 0.800  fractional weight for phases from previous global cycle
 -g 1.100  solvent gamma flipping factor
 -G 0.700  FOM threshold for initial tripeptides and chain building
 -i unset  no structure inversion
 -k   4.5  minimum height/sigma for revised heavy atom sites
 -l     2  space for  2000000 reflections
 -L     6  minimum number of residues/chain (if more than 3 chains)
 -m    10  cycles of density modification
 -o        omit residues from fragment to optimize CC
 -p unset  no phosphate search
 -q     7   alpha-helix search
 -r  3.00  map resolution (multiplies maximum indices)
 -s 0.600  solvent fraction
 -S  2.42  radius of sphere of influence
 -t  4.00  time for initial searches (-t6 or more if difficult)
 -u   500  MB allocatable memory for fragment optimization
 -U  0.00  abort if less than this % of fragment CA retained within 0.7A
 -v 0.000  density sharpening factor
 -w 0.200  weight for experimental phases after cycle 1
 -x unset  no phase and trace diagnostics 
 -y  1.80  highest resolution in Ang. for starting phases from model
 -z unset  do not optimize heavy atoms

 Space group: P 41 2 2   Allowed origin shift code:  9

    23385 Reflections read from file shelxe-input.hkl

    23385 Unique data, highest resolution =  2.200 Angstroms

 Anisotropic scaling: intensities multiplied by
  0.000624h^2 +0.000624k^2 -0.000100l^2 +0.000000kl +0.000000hl +0.000000hk

    72 Reflections with d > 2.400 and      0 in range  2.400 > d > 2.200 added
 Fourier grid =   128 x  128 x   29     0.000 <= z <= 0.125

    92 Point spherical net set up with radius 2.42A
    24 Extra Fourier layers will be generated


   210 Atoms read from PDB format file shelxe-input.pda

     0 NCS groups defined in .pda input file
  17 Residue blocks eliminated out of  42, increasing CC from  8.82 to 11.40%
  25 residues written to shelxe-input.pdo

 Overall CC between native Eobs and Ecalc (from fragment) = 13.93%
 Starting from fragment phases truncated to  1.800A
 <wt> = 0.295 for fragment phases
 <wt> = 0.295, Contrast = 0.196, Connect. = 0.617 for dens.mod. cycle 1
 <wt> = 0.296, Contrast = 0.244, Connect. = 0.628 for dens.mod. cycle 2
 <wt> = 0.297, Contrast = 0.265, Connect. = 0.642 for dens.mod. cycle 3
 <wt> = 0.299, Contrast = 0.279, Connect. = 0.652 for dens.mod. cycle 4
 <wt> = 0.300, Contrast = 0.289, Connect. = 0.659 for dens.mod. cycle 5
 <wt> = 0.300, Contrast = 0.298, Connect. = 0.665 for dens.mod. cycle 6
 <wt> = 0.300, Contrast = 0.305, Connect. = 0.671 for dens.mod. cycle 7
 <wt> = 0.300, Contrast = 0.311, Connect. = 0.676 for dens.mod. cycle 8
 <wt> = 0.300, Contrast = 0.316, Connect. = 0.679 for dens.mod. cycle 9
 <wt> = 0.300, Contrast = 0.320, Connect. = 0.682 for dens.mod. cycle 10

 NOGO map generated for regions about rotation axes (if any)

  1345 peaks > 0.5 sigma used to seed fragment search
 Space for about  350 unique residues taking solvent into account

   139 potential (overlapping) alpha-helix fragments found

   #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  
 
  56.0% starting CA preserved within 1.0A,  52.0% / 0.7A and  44.0% / 0.5A

 CC for partial structure against native data =  14.46 %

 ------------------------------------------------------------------------------

 Global autotracing cycle   2

 <wt> = 0.468 for fragment phases
 <wt> = 0.300, Contrast = 0.246, Connect. = 0.649 for dens.mod. cycle 1
 <wt> = 0.300, Contrast = 0.274, Connect. = 0.651 for dens.mod. cycle 2
 <wt> = 0.300, Contrast = 0.310, Connect. = 0.679 for dens.mod. cycle 3
 <wt> = 0.300, Contrast = 0.324, Connect. = 0.689 for dens.mod. cycle 4
 <wt> = 0.300, Contrast = 0.333, Connect. = 0.694 for dens.mod. cycle 5
 <wt> = 0.300, Contrast = 0.338, Connect. = 0.698 for dens.mod. cycle 6
 <wt> = 0.300, Contrast = 0.341, Connect. = 0.700 for dens.mod. cycle 7
 <wt> = 0.300, Contrast = 0.345, Connect. = 0.702 for dens.mod. cycle 8
 <wt> = 0.300, Contrast = 0.347, Connect. = 0.702 for dens.mod. cycle 9
 <wt> = 0.300, Contrast = 0.348, Connect. = 0.703 for dens.mod. cycle 10

 NOGO map generated for regions about rotation axes (if any)

  1331 peaks > 0.5 sigma used to seed fragment search
 Space for about  350 unique residues taking solvent into account

   206 potential (overlapping) alpha-helix fragments found

   #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  
 A:   13  24.554   2.434   0.819   0.545   0.833   0.500   CB 1.839
 B:   12  16.211   2.379   0.435   0.486   0.818   0.636   CB 1.775
 C:    9  12.906   1.744   0.682   0.769   0.875   0.500   CB 1.400
 D:   10   5.099   1.504   0.055   0.682   0.778   0.667   CB 1.165
 E:    8   4.370   1.356   0.103   0.746   0.857   0.571   CB 1.034
 F:    8   6.876   1.402   0.601   0.630   0.857   0.571   CB 1.196
 G:    7   8.645   1.489   0.359   0.869   1.000   0.833   CB 1.274
 H:   10   6.926   1.398   0.612   0.656   0.778   0.667   CB 1.139
 I:    9   6.779   1.317   0.545   0.806   0.875   0.875   CB 1.031
 J:    8   3.589   1.222   0.425   0.770   0.857   0.714   CB 0.716
 K:    6   4.672   1.162   0.427   0.900   1.000   1.000   O 0.900
 L:    7   7.308   1.237   0.956   0.910   1.000   0.833   CB 0.875
 M:    8   2.651   1.108   0.184   0.652   0.714   0.714   N 0.929
 N:    6   3.770   1.071   0.268   0.895   1.000   1.000   N 0.891
 O:    8   4.202   1.158   0.376   0.759   0.857   0.714   CB 0.924
 P:    6   3.146   0.949   0.332   0.895   1.000   1.000   C 0.798
 Q:    6   2.465   0.990   0.161   0.887   1.000   1.000   CB 0.691

   104 potential tripeptides employed

   #(CA)   CFOM   DenFit  NH...O  SecStr  Rama98  Rama80  
 R:   10   9.277   1.900   1.149   0.233   0.667   0.444   CB 1.550
 S:   22   6.717   1.791   0.549   0.149   0.524   0.238   CB 1.554
 T:    8   3.999   1.293   0.145   0.453   1.000   0.429   CB 1.071
       6   1.967   1.237   0.178   0.408   0.800   0.600   O 0.830 ?
 U:    8   3.333   1.624   0.353   0.284   0.857   0.429   CB 0.853
 V:   10   3.215   1.755   0.403   0.205   0.444   0.444   CB 1.390
 W:   14  12.759   2.349   0.371   0.378   0.769   0.385   CB 1.627
    16 residues pruned to eliminate duplicates
 X:    6   2.283   1.494  -0.022   0.347   0.800   0.400   CB 1.032
       6   1.950   1.143   0.400   0.379   0.800   0.600   CB 0.775 ?

 Using tripeptides from previous cycle as seeds

   141 residues left after pruning, divided into chains as follows:
 A:  12   B:   9   C:   8   D:   8   E:   7   F:  10   G:   9   H:   8   I:   7
 J:   8   K:   6   L:   6   M:   6   N:   6   O:  12   P:   8   Q:  11

  60.0% starting CA preserved within 1.0A,  56.0% / 0.7A and  32.0% / 0.5A

 CC for partial structure against native data =  17.52 %

  ------------------------------------------------------------------------------
 
 Global autotracing cycle  3

 Phases from autotracing cycle  2 used as input for final density modification

 <wt> = 0.538 for fragment phases
 <wt> = 0.300, Contrast = 0.239, Connect. = 0.665 for dens.mod. cycle 1
 <wt> = 0.300, Contrast = 0.268, Connect. = 0.668 for dens.mod. cycle 2
 <wt> = 0.300, Contrast = 0.314, Connect. = 0.708 for dens.mod. cycle 3
 <wt> = 0.300, Contrast = 0.331, Connect. = 0.720 for dens.mod. cycle 4
 <wt> = 0.300, Contrast = 0.341, Connect. = 0.722 for dens.mod. cycle 5
 <wt> = 0.300, Contrast = 0.344, Connect. = 0.723 for dens.mod. cycle 6
 <wt> = 0.300, Contrast = 0.349, Connect. = 0.724 for dens.mod. cycle 7
 <wt> = 0.300, Contrast = 0.350, Connect. = 0.724 for dens.mod. cycle 8
 <wt> = 0.300, Contrast = 0.353, Connect. = 0.723 for dens.mod. cycle 9
 <wt> = 0.300, Contrast = 0.355, Connect. = 0.723 for dens.mod. cycle 10

 Estimated mean FOM and mapCC as a function of resolution
 d    inf - 4.85 - 3.82 - 3.33 - 3.02 - 2.79 - 2.63 - 2.49 - 2.38 - 2.29 - 2.21
 <FOM>   0.526  0.502  0.439  0.406  0.358  0.286  0.310  0.339  0.315  0.265
 <mapCC> 0.942  0.863  0.765  0.753  0.674  0.547  0.593  0.637  0.480  0.389
 N        2347   2365   2321   2332   2448   2226   2452   2336   2259   2299

 Estimated mean FOM = 0.375   Pseudo-free CC = 43.82 %

 Best trace (cycle 3 with CC 25.64%) was saved as shelxe-input.pdb

 ==============================================================================

       CPU times required in seconds
       -----------------------------
         2.5 - Setup, data input and phasing
        25.3 - FFTs and peak-searches
        64.1 - Sphere of influence
         2.0 - Rest of density modification
      2515.9 - Alpha-helix search
      1218.5 - Tripeptide search
      1418.3 - Chain tracing
         0.0 - NCS analysis
       103.2 - B-value refinement for trace
         2.0 - Rest of tracing
         0.0 - Comparison with known structure


 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +  SHELXE finished at 04:31:58      Total time:      5351.91 secs  +
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 """
        file_contents = """TITLE  shelxe-input.pdb   Cycle  19   CC = 25.64%   183  residues in  14  chains
CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2
SCALE1      0.013637  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013637  0.000000        0.00000
SCALE3      0.000000  0.000000  0.006115        0.00000
ATOM      1  CA  ALA A   1      24.175  11.677 -28.178  1.00 20.00           C
ATOM      2  C   ALA A   1      22.865  12.325 -27.744  1.00 20.00           C
ATOM      3  O   ALA A   1      21.841  12.145 -28.400  1.00 20.00           O
ATOM      4  N   ALA A   2      22.893  13.079 -26.641  1.00 20.00           N
ATOM      5  CA  ALA A   2      21.697  13.745 -26.137  1.00 20.00           C
ATOM      6  CB  ALA A   2      21.427  15.034 -26.920  1.00 20.00           C
ATOM      7  C   ALA A   2      20.456  12.858 -26.179  1.00 20.00           C
ATOM      8  O   ALA A   2      19.348  13.357 -26.358  1.00 20.00           O
ATOM      9  N   ALA A   3      20.639  11.545 -26.010  1.00 20.00           N
ATOM     10  CA  ALA A   3      19.522  10.607 -26.031  1.00 20.00           C
ATOM    909  N   ALA N   6      12.748  28.352  30.277  1.00 20.00           N
ATOM    910  CA  ALA N   6      13.022  27.069  30.915  1.00 20.00           C
ATOM    911  CB  ALA N   6      14.129  26.296  30.191  1.00 20.00           C
ATOM    912  C   ALA N   6      13.407  27.321  32.369  1.00 20.00           C
ATOM    913  O   ALA N   6      12.957  26.603  33.259  1.00 20.00           O
ATOM    914  N   ALA N   7      14.237  28.339  32.613  1.00 20.00           N
ATOM    915  CA  ALA N   7      14.668  28.669  33.968  1.00 20.00           C
END 
"""

        fname = create_tempfile(content=file_contents)
        self.addCleanup(os.remove, fname)
        parser = ShelxeParser(fname=fname, stdout=stdout_contents)
        parser.parse()

        self.assertEqual("25.64", parser.cc)
        self.assertEqual("13.0", parser.acl)
        self.assertEqual("YES", parser.solution)
        self.assertEqual(3.0599999999999987, parser.average_cc_delta)
        self.assertEqual('13.93', parser.cc_eobs_ecalc)
        self.assertTupleEqual(('25.64', '13.0', '13.93', 3.0599999999999987, 'YES'), parser.summary)


if __name__ == '__main__':
    unittest.main()
