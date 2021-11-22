import os
import swamp
import unittest
from swamp.mr.searchmodel import SearchModel
from swamp.utils import remove, create_tempfile, compress, get_tempfile

TEST_FULLSIZE_PDB_SIDECHAINS = """HEADER    HYDROLASE/MEMBRANE PROTEIN              11-NOV-13   4NJN              
TITLE     CRYSTAL STRUCTURE OF E.COLI GLPG AT PH 4.5                            
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: RHOMBOID PROTEASE GLPG;                                    
COMPND   3 CHAIN: A;                                                            
COMPND   4 FRAGMENT: UNP RESIDUES 87-276;                                       
COMPND   5 SYNONYM: INTRAMEMBRANE SERINE PROTEASE;                              
COMPND   6 EC: 3.4.21.105;                                                      
COMPND   7 ENGINEERED: YES                                                      
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: ESCHERICHIA COLI;                               
SOURCE   3 ORGANISM_TAXID: 562;                                                 
SOURCE   4 GENE: GLPG, B3424, JW5687;                                           
SOURCE   5 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
SOURCE   6 EXPRESSION_SYSTEM_TAXID: 562                                         
KEYWDS    RHOMBOID PROTEASE, INTRAMEMBRANE PROTEOLYSIS, MEMBRANE PROTEIN,       
KEYWDS   2 HYDROLASE-MEMBRANE PROTEIN COMPLEX                                   
EXPDTA    X-RAY DIFFRACTION                                                     
AUTHOR    S.W.DICKEY,R.P.BAKER,S.CHO,S.URBAN                                    
REVDAT   2   01-JAN-14 4NJN    1       JRNL                                     
REVDAT   1   25-DEC-13 4NJN    0                                                
JRNL        AUTH   S.W.DICKEY,R.P.BAKER,S.CHO,S.URBAN                           
JRNL        TITL   PROTEOLYSIS INSIDE THE MEMBRANE IS A RATE-GOVERNED REACTION  
JRNL        TITL 2 NOT DRIVEN BY SUBSTRATE AFFINITY.                            
JRNL        REF    CELL(CAMBRIDGE,MASS.)         V. 155  1270 2013              
JRNL        REFN                   ISSN 0092-8674                               
JRNL        PMID   24315097                                                     
JRNL        DOI    10.1016/J.CELL.2013.10.053                                   
REMARK   2                                                                      
REMARK   2 RESOLUTION.    2.40 ANGSTROMS.                                       
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : REFMAC 5.7.0029                                      
REMARK   3   AUTHORS     : MURSHUDOV,VAGIN,DODSON                               
REMARK   3                                                                      
REMARK   3    REFINEMENT TARGET : MAXIMUM LIKELIHOOD                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.40                           
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 44.91                          
REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL                           
REMARK   3   COMPLETENESS FOR RANGE        (%) : 79.4                           
REMARK   3   NUMBER OF REFLECTIONS             : 9029                           
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.198                           
REMARK   3   R VALUE            (WORKING SET) : 0.196                           
REMARK   3   FREE R VALUE                     : 0.248                           
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.000                           
REMARK   3   FREE R VALUE TEST SET COUNT      : 478                             
REMARK   3                                                                      
REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  
REMARK   3   TOTAL NUMBER OF BINS USED           : 20                           
REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 2.40                         
REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.46                         
REMARK   3   REFLECTION IN BIN     (WORKING SET) : 389                          
REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 47.85                        
REMARK   3   BIN R VALUE           (WORKING SET) : 0.2590                       
REMARK   3   BIN FREE R VALUE SET COUNT          : 23                           
REMARK   3   BIN FREE R VALUE                    : 0.2920                       
REMARK   3                                                                      
REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    
REMARK   3   PROTEIN ATOMS            : 1451                                    
REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
REMARK   3   HETEROGEN ATOMS          : 0                                       
REMARK   3   SOLVENT ATOMS            : 50                                      
REMARK   3                                                                      
REMARK   3  B VALUES.                                                           
REMARK   3   FROM WILSON PLOT           (A**2) : NULL                           
REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 47.53                          
REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
REMARK   3    B11 (A**2) : -1.71000                                             
REMARK   3    B22 (A**2) : -1.71000                                             
REMARK   3    B33 (A**2) : 5.53000                                              
REMARK   3    B12 (A**2) : -1.71000                                             
REMARK   3    B13 (A**2) : -0.00000                                             
REMARK   3    B23 (A**2) : 0.00000                                              
REMARK   3                                                                      
REMARK   3  ESTIMATED OVERALL COORDINATE ERROR.                                 
REMARK   3   ESU BASED ON R VALUE                            (A): 0.364         
REMARK   3   ESU BASED ON FREE R VALUE                       (A): 0.266         
REMARK   3   ESU BASED ON MAXIMUM LIKELIHOOD                 (A): 0.155         
REMARK   3   ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): 6.591         
REMARK   3                                                                      
REMARK   3 CORRELATION COEFFICIENTS.                                            
REMARK   3   CORRELATION COEFFICIENT FO-FC      : 0.937                         
REMARK   3   CORRELATION COEFFICIENT FO-FC FREE : 0.908                         
REMARK   3                                                                      
REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT      
REMARK   3   BOND LENGTHS REFINED ATOMS        (A):  1501 ; 0.018 ; 0.019       
REMARK   3   BOND LENGTHS OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):  2043 ; 1.965 ; 1.924       
REMARK   3   BOND ANGLES OTHERS          (DEGREES):  NULL ;  NULL ;  NULL       
REMARK   3   TORSION ANGLES, PERIOD 1    (DEGREES):   181 ; 7.098 ; 5.000       
REMARK   3   TORSION ANGLES, PERIOD 2    (DEGREES):    59 ;30.230 ;22.203       
REMARK   3   TORSION ANGLES, PERIOD 3    (DEGREES):   233 ;18.954 ;15.000       
REMARK   3   TORSION ANGLES, PERIOD 4    (DEGREES):     6 ;16.340 ;15.000       
REMARK   3   CHIRAL-CENTER RESTRAINTS       (A**3):   219 ; 0.126 ; 0.200       
REMARK   3   GENERAL PLANES REFINED ATOMS      (A):  1127 ; 0.009 ; 0.021       
REMARK   3   GENERAL PLANES OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) REFINED ATOMS      (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW REFINED ATOMS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND REFINED ATOMS     (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND OTHERS            (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT      
REMARK   3   MAIN-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B REFINED ATOMS     (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B OTHER ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3 ANISOTROPIC THERMAL FACTOR RESTRAINTS.    COUNT   RMS   WEIGHT       
REMARK   3   RIGID-BOND RESTRAINTS          (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; FREE ATOMS         (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; BONDED ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  NCS RESTRAINTS STATISTICS                                           
REMARK   3   NUMBER OF DIFFERENT NCS GROUPS : NULL                              
REMARK   3                                                                      
REMARK   3  TLS DETAILS                                                         
REMARK   3   NUMBER OF TLS GROUPS  : NULL                                       
REMARK   3                                                                      
REMARK   3  BULK SOLVENT MODELLING.                                             
REMARK   3   METHOD USED : MASK                                                 
REMARK   3   PARAMETERS FOR MASK CALCULATION                                    
REMARK   3   VDW PROBE RADIUS   : 1.20                                          
REMARK   3   ION PROBE RADIUS   : 0.80                                          
REMARK   3   SHRINKAGE RADIUS   : 0.80                                          
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      
REMARK   4                                                                      
REMARK   4 4NJN COMPLIES WITH FORMAT V. 3.30, 13-JUL-11                             
CRYST1  110.680  110.680  127.710  90.00  90.00 120.00 H 3 2        18          
ORIGX1      1.000000  0.000000  0.000000        0.00000                         
ORIGX2      0.000000  1.000000  0.000000        0.00000                         
ORIGX3      0.000000  0.000000  1.000000        0.00000                         
SCALE1      0.009035  0.005216  0.000000        0.00000                         
SCALE2      0.000000  0.010433  0.000000        0.00000                         
SCALE3      0.000000  0.000000  0.007830        0.00000                         
ATOM      1  N   GLU A  91       7.010   5.663  26.232  1.00118.15           N  
ATOM      2  CA  GLU A  91       6.968   6.498  27.471  1.00124.45           C  
ATOM      3  C   GLU A  91       8.153   7.501  27.547  1.00114.42           C  
ATOM      4  O   GLU A  91       8.316   8.333  26.642  1.00 98.76           O  
ATOM      5  CB  GLU A  91       6.858   5.603  28.732  1.00130.03           C  
ATOM      6  CG  GLU A  91       5.452   5.386  29.300  1.00124.91           C  
ATOM      7  CD  GLU A  91       5.048   6.437  30.339  1.00127.68           C  
ATOM      8  OE1 GLU A  91       3.931   6.311  30.883  1.00121.99           O  
ATOM      9  OE2 GLU A  91       5.826   7.390  30.625  1.00114.91           O  
ATOM     10  N   ARG A  92       8.960   7.386  28.615  1.00112.39           N  
ATOM     11  CA  ARG A  92      10.003   8.358  29.069  1.00100.01           C  
ATOM     12  C   ARG A  92       9.446   9.729  29.467  1.00101.76           C  
ATOM     13  O   ARG A  92       8.724  10.352  28.692  1.00102.85           O  
ATOM     14  CB  ARG A  92      11.159   8.504  28.076  1.00 96.59           C  
ATOM     15  CG  ARG A  92      12.446   9.074  28.660  1.00 94.72           C  
ATOM     16  CD  ARG A  92      13.337   9.534  27.518  1.00 96.40           C  
ATOM     17  NE  ARG A  92      14.605  10.129  27.933  1.00 97.22           N  
ATOM     18  CZ  ARG A  92      15.505  10.638  27.088  1.00 99.40           C  
ATOM     19  NH1 ARG A  92      15.289  10.630  25.770  1.00 90.25           N  
ATOM     20  NH2 ARG A  92      16.631  11.157  27.561  1.00 92.33           N  
ATOM     21  N   ALA A  93       9.812  10.184  30.672  1.00 97.84           N  
ATOM     22  CA  ALA A  93       9.269  11.401  31.318  1.00 72.00           C  
ATOM     23  C   ALA A  93      10.202  12.617  31.152  1.00 66.40           C  
ATOM     24  O   ALA A  93      11.384  12.500  31.481  1.00 61.73           O  
ATOM     25  CB  ALA A  93       9.061  11.118  32.802  1.00 66.77           C  
ATOM     26  N   GLY A  94       9.663  13.770  30.709  1.00 59.55           N  
ATOM     27  CA  GLY A  94      10.430  15.030  30.447  1.00 49.91           C  
ATOM     28  C   GLY A  94      11.329  15.626  31.550  1.00 43.72           C  
ATOM     29  O   GLY A  94      11.204  15.252  32.750  1.00 46.03           O  
ATOM     30  N   PRO A  95      12.221  16.589  31.186  1.00 34.48           N  
ATOM     31  CA  PRO A  95      13.219  17.073  32.186  1.00 30.04           C  
ATOM     32  C   PRO A  95      12.651  17.690  33.487  1.00 30.74           C  
ATOM     33  O   PRO A  95      13.298  17.565  34.549  1.00 37.24           O  
ATOM     34  CB  PRO A  95      14.023  18.118  31.405  1.00 26.31           C  
ATOM     35  CG  PRO A  95      13.043  18.649  30.409  1.00 30.11           C  
ATOM     36  CD  PRO A  95      12.237  17.420  29.961  1.00 32.20           C  
ATOM     37  N   VAL A  96      11.531  18.397  33.422  1.00 27.66           N  
ATOM     38  CA  VAL A  96      10.984  19.038  34.645  1.00 30.64           C  
ATOM     39  C   VAL A  96      10.431  17.940  35.550  1.00 30.45           C  
ATOM     40  O   VAL A  96      10.788  17.877  36.749  1.00 28.11           O  
ATOM     41  CB  VAL A  96       9.876  20.107  34.375  1.00 29.89           C  
ATOM     42  CG1 VAL A  96       9.356  20.719  35.704  1.00 22.97           C  
ATOM     43  CG2 VAL A  96      10.405  21.169  33.383  1.00 27.98           C  
ATOM     44  N   THR A  97       9.608  17.059  34.961  1.00 28.09           N  
ATOM     45  CA  THR A  97       9.070  15.885  35.655  1.00 29.25           C  
ATOM     46  C   THR A  97      10.124  15.056  36.380  1.00 34.27           C  
ATOM     47  O   THR A  97       9.993  14.769  37.590  1.00 34.88           O  
ATOM     48  CB  THR A  97       8.355  14.996  34.661  1.00 35.14           C  
ATOM     49  OG1 THR A  97       7.368  15.794  33.959  1.00 38.68           O  
ATOM     50  CG2 THR A  97       7.693  13.838  35.383  1.00 25.85           C  
ATOM     51  N   TRP A  98      11.195  14.731  35.650  1.00 34.44           N  
ATOM     52  CA  TRP A  98      12.342  13.965  36.132  1.00 35.59           C  
ATOM     53  C   TRP A  98      13.210  14.670  37.198  1.00 38.80           C  
ATOM     54  O   TRP A  98      13.566  14.065  38.215  1.00 36.72           O  
ATOM     55  CB  TRP A  98      13.176  13.654  34.898  1.00 38.65           C  
ATOM     56  CG  TRP A  98      14.491  13.000  35.093  1.00 50.81           C  
ATOM     57  CD1 TRP A  98      15.703  13.463  34.653  1.00 58.06           C  
ATOM     58  CD2 TRP A  98      14.746  11.733  35.702  1.00 55.64           C  
ATOM     59  NE1 TRP A  98      16.691  12.570  34.957  1.00 58.98           N  
ATOM     60  CE2 TRP A  98      16.136  11.503  35.614  1.00 62.29           C  
ATOM     61  CE3 TRP A  98      13.937  10.769  36.327  1.00 60.49           C  
ATOM     62  CZ2 TRP A  98      16.745  10.347  36.146  1.00 62.49           C  
ATOM     63  CZ3 TRP A  98      14.538   9.623  36.842  1.00 63.53           C  
ATOM     64  CH2 TRP A  98      15.932   9.421  36.746  1.00 57.46           C  
ATOM     65  N   VAL A  99      13.563  15.943  36.972  1.00 38.35           N  
ATOM     66  CA  VAL A  99      14.410  16.636  37.926  1.00 37.34           C  
ATOM     67  C   VAL A  99      13.748  16.752  39.279  1.00 34.52           C  
ATOM     68  O   VAL A  99      14.391  16.590  40.283  1.00 37.41           O  
ATOM     69  CB  VAL A  99      14.910  17.998  37.413  1.00 41.19           C  
ATOM     70  CG1 VAL A  99      15.572  18.796  38.542  1.00 38.55           C  
ATOM     71  CG2 VAL A  99      15.916  17.740  36.293  1.00 40.93           C  
ATOM     72  N   MET A 100      12.448  16.974  39.275  1.00 34.10           N  
ATOM     73  CA  MET A 100      11.685  17.117  40.485  1.00 36.21           C  
ATOM     74  C   MET A 100      11.571  15.802  41.309  1.00 38.30           C  
ATOM     75  O   MET A 100      11.427  15.843  42.564  1.00 33.46           O  
ATOM     76  CB  MET A 100      10.300  17.605  40.101  1.00 36.03           C  
ATOM     77  CG  MET A 100       9.332  17.538  41.264  1.00 52.58           C  
ATOM     78  SD  MET A 100       9.615  18.768  42.569  1.00 63.73           S  
ATOM     79  CE  MET A 100       8.818  20.163  41.773  1.00 60.77           C  
MASTER      444    0    0   11    0    0    0    6 1501    1    0   17          
END                                                                             
"""

TEST_FULLSIZE_PDB_POLYALA = """HEADER    HYDROLASE/MEMBRANE PROTEIN              11-NOV-13   4NJN              
TITLE     CRYSTAL STRUCTURE OF E.COLI GLPG AT PH 4.5                            
KEYWDS    RHOMBOID PROTEASE, INTRAMEMBRANE PROTEOLYSIS, MEMBRANE PROTEIN,       
KEYWDS   2 HYDROLASE-MEMBRANE PROTEIN COMPLEX                                   
EXPDTA    X-RAY DIFFRACTION                                                     
REMARK   2                                                                      
REMARK   2 RESOLUTION.    2.40 ANGSTROMS.                                       
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : REFMAC 5.7.0029                                      
REMARK   3   AUTHORS     : MURSHUDOV,VAGIN,DODSON                               
REMARK   3                                                                      
REMARK   3    REFINEMENT TARGET : MAXIMUM LIKELIHOOD                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.40                           
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 44.91                          
REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL                           
REMARK   3   COMPLETENESS FOR RANGE        (%) : 79.4                           
REMARK   3   NUMBER OF REFLECTIONS             : 9029                           
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.198                           
REMARK   3   R VALUE            (WORKING SET) : 0.196                           
REMARK   3   FREE R VALUE                     : 0.248                           
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.000                           
REMARK   3   FREE R VALUE TEST SET COUNT      : 478                             
REMARK   3                                                                      
REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  
REMARK   3   TOTAL NUMBER OF BINS USED           : 20                           
REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 2.40                         
REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.46                         
REMARK   3   REFLECTION IN BIN     (WORKING SET) : 389                          
REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 47.85                        
REMARK   3   BIN R VALUE           (WORKING SET) : 0.2590                       
REMARK   3   BIN FREE R VALUE SET COUNT          : 23                           
REMARK   3   BIN FREE R VALUE                    : 0.2920                       
REMARK   3                                                                      
REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    
REMARK   3   PROTEIN ATOMS            : 1451                                    
REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
REMARK   3   HETEROGEN ATOMS          : 0                                       
REMARK   3   SOLVENT ATOMS            : 50                                      
REMARK   3                                                                      
REMARK   3  B VALUES.                                                           
REMARK   3   FROM WILSON PLOT           (A**2) : NULL                           
REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 47.53                          
REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
REMARK   3    B11 (A**2) : -1.71000                                             
REMARK   3    B22 (A**2) : -1.71000                                             
REMARK   3    B33 (A**2) : 5.53000                                              
REMARK   3    B12 (A**2) : -1.71000                                             
REMARK   3    B13 (A**2) : -0.00000                                             
REMARK   3    B23 (A**2) : 0.00000                                              
REMARK   3                                                                      
REMARK   3  ESTIMATED OVERALL COORDINATE ERROR.                                 
REMARK   3   ESU BASED ON R VALUE                            (A): 0.364         
REMARK   3   ESU BASED ON FREE R VALUE                       (A): 0.266         
REMARK   3   ESU BASED ON MAXIMUM LIKELIHOOD                 (A): 0.155         
REMARK   3   ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): 6.591         
REMARK   3                                                                      
REMARK   3 CORRELATION COEFFICIENTS.                                            
REMARK   3   CORRELATION COEFFICIENT FO-FC      : 0.937                         
REMARK   3   CORRELATION COEFFICIENT FO-FC FREE : 0.908                         
REMARK   3                                                                      
REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT      
REMARK   3   BOND LENGTHS REFINED ATOMS        (A):  1501 ; 0.018 ; 0.019       
REMARK   3   BOND LENGTHS OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):  2043 ; 1.965 ; 1.924       
REMARK   3   BOND ANGLES OTHERS          (DEGREES):  NULL ;  NULL ;  NULL       
REMARK   3   TORSION ANGLES, PERIOD 1    (DEGREES):   181 ; 7.098 ; 5.000       
REMARK   3   TORSION ANGLES, PERIOD 2    (DEGREES):    59 ;30.230 ;22.203       
REMARK   3   TORSION ANGLES, PERIOD 3    (DEGREES):   233 ;18.954 ;15.000       
REMARK   3   TORSION ANGLES, PERIOD 4    (DEGREES):     6 ;16.340 ;15.000       
REMARK   3   CHIRAL-CENTER RESTRAINTS       (A**3):   219 ; 0.126 ; 0.200       
REMARK   3   GENERAL PLANES REFINED ATOMS      (A):  1127 ; 0.009 ; 0.021       
REMARK   3   GENERAL PLANES OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) REFINED ATOMS      (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW REFINED ATOMS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND REFINED ATOMS     (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND OTHERS            (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT      
REMARK   3   MAIN-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B REFINED ATOMS     (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B OTHER ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3 ANISOTROPIC THERMAL FACTOR RESTRAINTS.    COUNT   RMS   WEIGHT       
REMARK   3   RIGID-BOND RESTRAINTS          (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; FREE ATOMS         (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; BONDED ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  NCS RESTRAINTS STATISTICS                                           
REMARK   3   NUMBER OF DIFFERENT NCS GROUPS : NULL                              
REMARK   3                                                                      
REMARK   3  TLS DETAILS                                                         
REMARK   3   NUMBER OF TLS GROUPS  : NULL                                       
REMARK   3                                                                      
REMARK   3  BULK SOLVENT MODELLING.                                             
REMARK   3   METHOD USED : MASK                                                 
REMARK   3   PARAMETERS FOR MASK CALCULATION                                    
REMARK   3   VDW PROBE RADIUS   : 1.20                                          
REMARK   3   ION PROBE RADIUS   : 0.80                                          
REMARK   3   SHRINKAGE RADIUS   : 0.80                                          
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      
REMARK   4                                                                      
REMARK   4 4NJN COMPLIES WITH FORMAT V. 3.30, 13-JUL-11                             
CRYST1  110.680  110.680  127.710  90.00  90.00 120.00 H 3 2        18          
ATOM      1  N   ALA A   1       7.010   5.663  26.232  1.00118.15           N  
ATOM      2  CA  ALA A   1       6.968   6.498  27.471  1.00124.45           C  
ATOM      3  C   ALA A   1       8.153   7.501  27.547  1.00114.42           C  
ATOM      4  O   ALA A   1       8.316   8.333  26.642  1.00 98.76           O  
ATOM      5  CB  ALA A   1       6.858   5.603  28.732  1.00130.03           C  
ATOM      6  N   ALA A   2       8.960   7.386  28.615  1.00112.39           N  
ATOM      7  CA  ALA A   2      10.003   8.358  29.069  1.00100.01           C  
ATOM      8  C   ALA A   2       9.446   9.729  29.467  1.00101.76           C  
ATOM      9  O   ALA A   2       8.724  10.352  28.692  1.00102.85           O  
ATOM     10  CB  ALA A   2      11.159   8.504  28.076  1.00 96.59           C  
ATOM     11  N   ALA A   3       9.812  10.184  30.672  1.00 97.84           N  
ATOM     12  CA  ALA A   3       9.269  11.401  31.318  1.00 72.00           C  
ATOM     13  C   ALA A   3      10.202  12.617  31.152  1.00 66.40           C  
ATOM     14  O   ALA A   3      11.384  12.500  31.481  1.00 61.73           O  
ATOM     15  CB  ALA A   3       9.061  11.118  32.802  1.00 66.77           C  
ATOM     16  N   ALA A   4       9.663  13.770  30.709  1.00 59.55           N  
ATOM     17  CA  ALA A   4      10.430  15.030  30.447  1.00 49.91           C  
ATOM     18  C   ALA A   4      11.329  15.626  31.550  1.00 43.72           C  
ATOM     19  O   ALA A   4      11.204  15.252  32.750  1.00 46.03           O  
ATOM     20  N   ALA A   5      12.221  16.589  31.186  1.00 34.48           N  
ATOM     21  CA  ALA A   5      13.219  17.073  32.186  1.00 30.04           C  
ATOM     22  C   ALA A   5      12.651  17.690  33.487  1.00 30.74           C  
ATOM     23  O   ALA A   5      13.298  17.565  34.549  1.00 37.24           O  
ATOM     24  CB  ALA A   5      14.023  18.118  31.405  1.00 26.31           C  
ATOM     25  N   ALA A   6      11.531  18.397  33.422  1.00 27.66           N  
ATOM     26  CA  ALA A   6      10.984  19.038  34.645  1.00 30.64           C  
ATOM     27  C   ALA A   6      10.431  17.940  35.550  1.00 30.45           C  
ATOM     28  O   ALA A   6      10.788  17.877  36.749  1.00 28.11           O  
ATOM     29  CB  ALA A   6       9.876  20.107  34.375  1.00 29.89           C  
ATOM     30  N   ALA A   7       9.608  17.059  34.961  1.00 28.09           N  
ATOM     31  CA  ALA A   7       9.070  15.885  35.655  1.00 29.25           C  
ATOM     32  C   ALA A   7      10.124  15.056  36.380  1.00 34.27           C  
ATOM     33  O   ALA A   7       9.993  14.769  37.590  1.00 34.88           O  
ATOM     34  CB  ALA A   7       8.355  14.996  34.661  1.00 35.14           C  
ATOM     35  N   ALA A   8      11.195  14.731  35.650  1.00 34.44           N  
ATOM     36  CA  ALA A   8      12.342  13.965  36.132  1.00 35.59           C  
ATOM     37  C   ALA A   8      13.210  14.670  37.198  1.00 38.80           C  
ATOM     38  O   ALA A   8      13.566  14.065  38.215  1.00 36.72           O  
ATOM     39  CB  ALA A   8      13.176  13.654  34.898  1.00 38.65           C  
ATOM     40  N   ALA A   9      13.563  15.943  36.972  1.00 38.35           N  
ATOM     41  CA  ALA A   9      14.410  16.636  37.926  1.00 37.34           C  
ATOM     42  C   ALA A   9      13.748  16.752  39.279  1.00 34.52           C  
ATOM     43  O   ALA A   9      14.391  16.590  40.283  1.00 37.41           O  
ATOM     44  CB  ALA A   9      14.910  17.998  37.413  1.00 41.19           C  
ATOM     45  N   ALA A  10      12.448  16.974  39.275  1.00 34.10           N  
ATOM     46  CA  ALA A  10      11.685  17.117  40.485  1.00 36.21           C  
ATOM     47  C   ALA A  10      11.571  15.802  41.309  1.00 38.30           C  
ATOM     48  O   ALA A  10      11.427  15.843  42.564  1.00 33.46           O  
ATOM     49  CB  ALA A  10      10.300  17.605  40.101  1.00 36.03           C  
END  """

ENSEMBLE_SIDECHAINS = """NUMMDL    2                                                                     
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1                      
MODEL        1                                                                  
ATOM      1  N   VAL A   1       6.808  61.263  13.432  1.00 72.61           N  
ANISOU    1  N   VAL A   1     9095   7963  10530    149   1518  -1669       N  
ATOM      2  CA  VAL A   1       7.802  60.901  14.433  1.00 75.34           C  
ANISOU    2  CA  VAL A   1     9317   8124  11185    233   1511  -1496       C  
ATOM      3  C   VAL A   1       8.710  59.751  13.957  1.00 83.54           C  
ANISOU    3  C   VAL A   1    10223   8997  12519    336   1668  -1701       C  
ATOM      4  O   VAL A   1       9.847  59.594  14.417  1.00 89.53           O  
ANISOU    4  O   VAL A   1    10845   9666  13503    427   1706  -1596       O  
ATOM      5  CB  VAL A   1       7.145  60.691  15.789  1.00 75.55           C  
ANISOU    5  CB  VAL A   1     9364   7997  11342    213   1344  -1273       C  
ATOM      6  CG1 VAL A   1       8.209  60.606  16.857  1.00 75.82           C  
ANISOU    6  CG1 VAL A   1     9281   7911  11614    291   1304  -1042       C  
ATOM      7  CG2 VAL A   1       6.221  61.872  16.098  1.00 67.42           C  
ANISOU    7  CG2 VAL A   1     8458   7155  10003    125   1217  -1130       C  
ATOM      8  N   ARG A   2       5.296  61.322  10.847  1.00 71.18           N  
ANISOU    8  N   ARG A   2     9096   8147   9801     20   1606  -2181       N  
ATOM      9  CA  ARG A   2       6.152  62.390  11.390  1.00 71.60           C  
ANISOU    9  CA  ARG A   2     9122   8276   9806     47   1598  -1902       C  
ATOM     10  C   ARG A   2       7.216  61.908  12.350  1.00 76.14           C  
ANISOU   10  C   ARG A   2     9562   8635  10730    135   1627  -1774       C  
ATOM     11  O   ARG A   2       8.425  62.127  12.148  1.00 78.15           O  
ANISOU   11  O   ARG A   2     9721   8935  11036    188   1743  -1750       O  
ATOM     12  CB  ARG A   2       5.322  63.389  12.190  1.00 81.86           C  
ANISOU   12  CB  ARG A   2    10508   9630  10963     -5   1414  -1647       C  
ATOM     13  CG  ARG A   2       4.263  64.128  11.397  1.00 81.49           C  
ANISOU   13  CG  ARG A   2    10585   9807  10569    -78   1344  -1694       C  
ATOM     14  CD  ARG A   2       4.883  64.688  10.150  1.00 80.59           C  
ANISOU   14  CD  ARG A   2    10495   9923  10201    -90   1470  -1776       C  
ATOM     15  NE  ARG A   2       4.395  66.021   9.919  1.00 78.92           N  
ANISOU   15  NE  ARG A   2    10386   9906   9693   -138   1372  -1604       N  
ATOM     16  CZ  ARG A   2       4.855  66.799   8.954  1.00 80.61           C  
ANISOU   16  CZ  ARG A   2    10641  10337   9649   -166   1452  -1579       C  
ATOM     17  NH1 ARG A   2       5.864  66.376   8.141  1.00 69.58           N  
ANISOU   17  NH1 ARG A   2     9180   9018   8239   -152   1651  -1730       N  
ATOM     18  NH2 ARG A   2       4.297  68.008   8.828  1.00 76.46           N  
ANISOU   18  NH2 ARG A   2    10215   9946   8888   -205   1335  -1395       N  
ATOM     19  N   ASP A   3       4.070  58.691  10.381  1.00 78.09           N  
ANISOU   19  N   ASP A   3     9947   8649  11074    -14   1625  -2713       N  
ATOM     20  CA  ASP A   3       4.632  59.476   9.274  1.00 77.89           C  
ANISOU   20  CA  ASP A   3     9952   8925  10717     -7   1737  -2805       C  
ATOM     21  C   ASP A   3       5.589  60.662   9.714  1.00 79.72           C  
ANISOU   21  C   ASP A   3    10158   9278  10852     34   1755  -2495       C  
ATOM     22  O   ASP A   3       6.580  60.956   9.019  1.00 79.47           O  
ANISOU   22  O   ASP A   3    10082   9389  10722     72   1910  -2556       O  
ATOM     23  CB  ASP A   3       3.483  59.989   8.383  1.00 78.46           C  
ANISOU   23  CB  ASP A   3    10144   9263  10404   -112   1654  -2920       C  
ATOM     24  CG  ASP A   3       2.890  58.904   7.491  1.00 87.78           C  
ANISOU   24  CG  ASP A   3    11337  10412  11602   -159   1689  -3319       C  
ATOM     25  OD1 ASP A   3       3.237  57.713   7.684  1.00 85.97           O  
ANISOU   25  OD1 ASP A   3    11032   9904  11728   -113   1765  -3497       O  
ATOM     26  OD2 ASP A   3       2.084  59.226   6.566  1.00 93.19           O1-
ANISOU   26  OD2 ASP A   3    12106  11348  11952   -240   1634  -3466       O1-
ATOM     27  N   PHE A   4       3.150  58.127  13.054  1.00 73.86           N  
ANISOU   27  N   PHE A   4     9370   7685  11008    -44   1372  -2230       N  
ATOM     28  CA  PHE A   4       4.164  57.312  12.379  1.00 79.41           C  
ANISOU   28  CA  PHE A   4     9998   8267  11906     42   1533  -2465       C  
ATOM     29  C   PHE A   4       4.849  58.044  11.235  1.00 79.34           C  
ANISOU   29  C   PHE A   4    10005   8539  11602     70   1658  -2607       C  
ATOM     30  O   PHE A   4       6.064  58.036  11.132  1.00 82.59           O  
ANISOU   30  O   PHE A   4    10329   8935  12116    166   1782  -2605       O  
ATOM     31  CB  PHE A   4       3.529  56.054  11.792  1.00 90.18           C  
ANISOU   31  CB  PHE A   4    11362   9448  13451     -3   1563  -2793       C  
ATOM     32  CG  PHE A   4       3.500  54.890  12.734  1.00 94.10           C  
ANISOU   32  CG  PHE A   4    11791   9558  14403     16   1528  -2712       C  
ATOM     33  CD1 PHE A   4       2.708  54.925  13.874  1.00 95.62           C  
ANISOU   33  CD1 PHE A   4    12002   9662  14668    -55   1384  -2429       C  
ATOM     34  CD2 PHE A   4       4.274  53.774  12.487  1.00 95.93           C  
ANISOU   34  CD2 PHE A   4    11939   9517  14992    111   1645  -2909       C  
ATOM     35  CE1 PHE A   4       2.687  53.872  14.754  1.00 93.35           C  
ANISOU   35  CE1 PHE A   4    11656   9027  14784    -47   1352  -2312       C  
ATOM     36  CE2 PHE A   4       4.258  52.709  13.358  1.00103.24           C  
ANISOU   36  CE2 PHE A   4    12805  10062  16356    132   1602  -2801       C  
ATOM     37  CZ  PHE A   4       3.459  52.759  14.492  1.00105.77           C  
ANISOU   37  CZ  PHE A   4    13153  10306  16728     45   1453  -2487       C  
ATOM     38  N   MET A   5       1.211  60.152  13.195  1.00 76.46           N  
ANISOU   38  N   MET A   5     9861   8451  10738   -188   1133  -1999       N  
ATOM     39  CA  MET A   5       2.272  60.073  14.171  1.00 72.86           C  
ANISOU   39  CA  MET A   5     9339   7838  10506   -107   1163  -1788       C  
ATOM     40  C   MET A   5       3.430  59.310  13.580  1.00 72.45           C  
ANISOU   40  C   MET A   5     9214   7668  10645    -28   1316  -1965       C  
ATOM     41  O   MET A   5       4.556  59.785  13.621  1.00 73.23           O  
ANISOU   41  O   MET A   5     9266   7815  10740     44   1382  -1871       O  
ATOM     42  CB  MET A   5       1.830  59.378  15.430  1.00 75.99           C  
ANISOU   42  CB  MET A   5     9696   8009  11167   -128   1083  -1618       C  
ATOM     43  CG  MET A   5       2.992  59.211  16.399  1.00 78.51           C  
ANISOU   43  CG  MET A   5     9940   8173  11717    -38   1099  -1400       C  
ATOM     44  SD  MET A   5       2.401  58.355  17.857  1.00 94.49           S  
ANISOU   44  SD  MET A   5    11929   9956  14017    -76   1002  -1172       S  
ATOM     45  CE  MET A   5       2.005  56.777  17.050  1.00 79.80           C  
ANISOU   45  CE  MET A   5    10036   7836  12447   -112   1080  -1481       C  
ATOM     46  N   ALA A   6      -0.524  59.587  11.161  1.00 66.38           N  
ANISOU   46  N   ALA A   6     8660   7395   9166   -337   1099  -2523       N  
ATOM     47  CA  ALA A   6       0.233  60.855  11.098  1.00 63.60           C  
ANISOU   47  CA  ALA A   6     8346   7236   8581   -274   1122  -2331       C  
ATOM     48  C   ALA A   6       1.411  60.765  12.032  1.00 67.49           C  
ANISOU   48  C   ALA A   6     8775   7548   9320   -188   1188  -2138       C  
ATOM     49  O   ALA A   6       2.505  61.188  11.690  1.00 69.06           O  
ANISOU   49  O   ALA A   6     8959   7820   9461   -127   1289  -2114       O  
ATOM     50  CB  ALA A   6      -0.614  62.070  11.372  1.00 56.93           C  
ANISOU   50  CB  ALA A   6     7562   6584   7483   -304    980  -2128       C  
ATOM     51  N   ILE A   7      -1.601  57.282  11.899  1.00 61.88           N  
ANISOU   51  N   ILE A   7     7992   6336   9182   -453   1058  -2717       N  
ATOM     52  CA  ILE A   7      -0.766  57.188  10.752  1.00 67.84           C  
ANISOU   52  CA  ILE A   7     8764   7177   9832   -397   1183  -2972       C  
ATOM     53  C   ILE A   7       0.013  58.494  10.624  1.00 70.11           C  
ANISOU   53  C   ILE A   7     9092   7710   9835   -319   1216  -2784       C  
ATOM     54  O   ILE A   7       1.136  58.503  10.080  1.00 69.35           O  
ANISOU   54  O   ILE A   7     8979   7645   9722   -241   1355  -2874       O  
ATOM     55  CB  ILE A   7      -1.441  56.665   9.450  1.00 73.01           C  
ANISOU   55  CB  ILE A   7     9448   7939  10350   -483   1185  -3372       C  
ATOM     56  CG1 ILE A   7      -1.949  57.797   8.612  1.00 74.68           C  
ANISOU   56  CG1 ILE A   7     9736   8541  10096   -517   1115  -3381       C  
ATOM     57  CG2 ILE A   7      -2.495  55.560   9.671  1.00 69.96           C  
ANISOU   57  CG2 ILE A   7     9022   7330  10227   -602   1106  -3522       C  
ATOM     58  CD1 ILE A   7      -0.901  58.278   7.651  1.00 81.84           C  
ANISOU   58  CD1 ILE A   7    10682   9651  10760   -443   1252  -3481       C  
ATOM     59  N   ILE A   8      -2.856  58.098  14.278  1.00 60.88           N  
ANISOU   59  N   ILE A   8     7836   6188   9105   -513    856  -2123       N  
ATOM     60  CA  ILE A   8      -1.791  57.132  14.351  1.00 59.72           C  
ANISOU   60  CA  ILE A   8     7651   5777   9263   -452    963  -2180       C  
ATOM     61  C   ILE A   8      -0.982  57.170  13.066  1.00 60.98           C  
ANISOU   61  C   ILE A   8     7835   6031   9301   -391   1071  -2449       C  
ATOM     62  O   ILE A   8       0.223  57.129  13.131  1.00 62.29           O  
ANISOU   62  O   ILE A   8     7974   6123   9569   -287   1166  -2409       O  
ATOM     63  CB  ILE A   8      -2.352  55.726  14.631  1.00 62.50           C  
ANISOU   63  CB  ILE A   8     7947   5828   9972   -539    960  -2276       C  
ATOM     64  CG1 ILE A   8      -2.751  55.640  16.091  1.00 58.98           C  
ANISOU   64  CG1 ILE A   8     7465   5266   9679   -575    893  -1936       C  
ATOM     65  CG2 ILE A   8      -1.344  54.630  14.376  1.00 63.71           C  
ANISOU   65  CG2 ILE A   8     8060   5689  10455   -468   1072  -2429       C  
ATOM     66  CD1 ILE A   8      -1.579  55.717  17.008  1.00 56.55           C  
ANISOU   66  CD1 ILE A   8     7137   4840   9506   -459    926  -1672       C  
ATOM     67  N   SER A   9      -4.846  59.990  13.395  1.00 58.40           N  
ANISOU   67  N   SER A   9     7589   6436   8161   -595    652  -2149       N  
ATOM     68  CA  SER A   9      -3.712  60.360  14.217  1.00 60.43           C  
ANISOU   68  CA  SER A   9     7859   6604   8495   -498    713  -1918       C  
ATOM     69  C   SER A   9      -2.570  59.381  14.182  1.00 58.45           C  
ANISOU   69  C   SER A   9     7580   6117   8509   -453    834  -1997       C  
ATOM     70  O   SER A   9      -1.430  59.784  14.089  1.00 59.36           O  
ANISOU   70  O   SER A   9     7710   6250   8590   -363    905  -1939       O  
ATOM     71  CB  SER A   9      -4.092  60.532  15.675  1.00 58.59           C  
ANISOU   71  CB  SER A   9     7590   6300   8370   -505    657  -1640       C  
ATOM     72  OG  SER A   9      -4.922  61.633  15.773  1.00 56.29           O  
ANISOU   72  OG  SER A   9     7321   6232   7833   -506    562  -1552       O  
ATOM     73  N   CYS A  10      -6.665  58.234  12.008  1.00 64.10           N  
ANISOU   73  N   CYS A  10     8232   7113   9010   -830    573  -2697       N  
ATOM     74  CA  CYS A  10      -5.894  59.264  11.295  1.00 66.28           C  
ANISOU   74  CA  CYS A  10     8594   7613   8975   -726    604  -2683       C  
ATOM     75  C   CYS A  10      -4.683  59.763  12.098  1.00 63.40           C  
ANISOU   75  C   CYS A  10     8255   7161   8671   -608    688  -2424       C  
ATOM     76  O   CYS A  10      -3.609  59.922  11.539  1.00 68.98           O  
ANISOU   76  O   CYS A  10     9001   7896   9309   -536    786  -2482       O  
ATOM     77  CB  CYS A  10      -6.774  60.431  10.945  1.00 62.41           C  
ANISOU   77  CB  CYS A  10     8127   7431   8153   -729    472  -2610       C  
ATOM     78  SG  CYS A  10      -8.209  59.959  10.007  1.00 70.01           S  
ANISOU   78  SG  CYS A  10     9038   8542   9019   -865    340  -2895       S  
ATOM     79  N   GLY A  11      -7.523  56.981  14.399  1.00 63.45           N  
ANISOU   79  N   GLY A  11     8001   6600   9506   -957    560  -2345       N  
ATOM     80  CA  GLY A  11      -6.926  56.234  13.306  1.00 64.48           C  
ANISOU   80  CA  GLY A  11     8161   6628   9707   -955    625  -2661       C  
ATOM     81  C   GLY A  11      -6.072  57.139  12.464  1.00 63.66           C  
ANISOU   81  C   GLY A  11     8144   6743   9299   -835    664  -2722       C  
ATOM     82  O   GLY A  11      -4.886  56.872  12.245  1.00 60.49           O  
ANISOU   82  O   GLY A  11     7771   6229   8982   -747    777  -2775       O  
ATOM     83  N   LEU A  12      -8.666  58.912  15.975  1.00 61.98           N  
ANISOU   83  N   LEU A  12     7779   6763   9007   -914    436  -1866       N  
ATOM     84  CA  LEU A  12      -7.453  58.270  16.446  1.00 60.60           C  
ANISOU   84  CA  LEU A  12     7634   6335   9053   -862    535  -1781       C  
ATOM     85  C   LEU A  12      -6.762  57.549  15.313  1.00 64.32           C  
ANISOU   85  C   LEU A  12     8139   6692   9606   -853    596  -2064       C  
ATOM     86  O   LEU A  12      -5.548  57.496  15.275  1.00 65.65           O  
ANISOU   86  O   LEU A  12     8348   6762   9831   -754    676  -2042       O  
ATOM     87  CB  LEU A  12      -7.713  57.288  17.567  1.00 60.13           C  
ANISOU   87  CB  LEU A  12     7504   6036   9306   -947    562  -1617       C  
ATOM     88  CG  LEU A  12      -6.496  56.752  18.325  1.00 60.68           C  
ANISOU   88  CG  LEU A  12     7594   5857   9604   -874    637  -1439       C  
ATOM     89  CD1 LEU A  12      -5.578  57.873  18.811  1.00 52.44           C  
ANISOU   89  CD1 LEU A  12     6610   4952   8360   -733    639  -1242       C  
ATOM     90  CD2 LEU A  12      -6.989  55.847  19.476  1.00 56.47           C  
ANISOU   90  CD2 LEU A  12     6989   5121   9344   -980    644  -1229       C  
ATOM     91  N   ALA A  13     -10.929  59.546  14.435  1.00 68.55           N  
ANISOU   91  N   ALA A  13     8519   8026   9501  -1036    226  -2184       N  
ATOM     92  CA  ALA A  13      -9.869  60.510  14.575  1.00 63.92           C  
ANISOU   92  CA  ALA A  13     8040   7496   8748   -885    267  -2029       C  
ATOM     93  C   ALA A  13      -8.582  59.839  15.040  1.00 63.57           C  
ANISOU   93  C   ALA A  13     8042   7192   8919   -846    392  -1965       C  
ATOM     94  O   ALA A  13      -7.516  60.131  14.527  1.00 72.56           O  
ANISOU   94  O   ALA A  13     9262   8339   9969   -758    447  -2000       O  
ATOM     95  CB  ALA A  13     -10.261  61.598  15.526  1.00 62.67           C  
ANISOU   95  CB  ALA A  13     7865   7469   8476   -815    225  -1776       C  
ATOM     96  N   TYR A  14     -12.378  57.236  14.902  1.00 66.70           N  
ANISOU   96  N   TYR A  14     8073   7456   9813  -1356    227  -2355       N  
ATOM     97  CA  TYR A  14     -11.949  57.558  13.541  1.00 67.37           C  
ANISOU   97  CA  TYR A  14     8246   7680   9670  -1296    192  -2604       C  
ATOM     98  C   TYR A  14     -10.808  58.532  13.580  1.00 65.65           C  
ANISOU   98  C   TYR A  14     8152   7540   9251  -1118    245  -2455       C  
ATOM     99  O   TYR A  14      -9.829  58.360  12.875  1.00 64.93           O  
ANISOU   99  O   TYR A  14     8148   7399   9123  -1062    310  -2590       O  
ATOM    100  CB  TYR A  14     -13.052  58.190  12.723  1.00 66.46           C  
ANISOU  100  CB  TYR A  14     8075   7875   9300  -1332     45  -2738       C  
ATOM    101  CG  TYR A  14     -12.714  58.385  11.254  1.00 64.29           C  
ANISOU  101  CG  TYR A  14     7888   7765   8773  -1300      0  -3004       C  
ATOM    102  CD1 TYR A  14     -12.713  57.314  10.356  1.00 65.34           C  
ANISOU  102  CD1 TYR A  14     8019   7800   9006  -1411      2  -3344       C  
ATOM    103  CD2 TYR A  14     -12.410  59.636  10.761  1.00 64.98           C  
ANISOU  103  CD2 TYR A  14     8062   8108   8517  -1166    -44  -2919       C  
ATOM    104  CE1 TYR A  14     -12.434  57.498   9.009  1.00 66.37           C  
ANISOU  104  CE1 TYR A  14     8229   8124   8862  -1388    -33  -3598       C  
ATOM    105  CE2 TYR A  14     -12.126  59.831   9.409  1.00 65.91           C  
ANISOU  105  CE2 TYR A  14     8261   8410   8369  -1148    -81  -3136       C  
ATOM    106  CZ  TYR A  14     -12.158  58.762   8.538  1.00 67.58           C  
ANISOU  106  CZ  TYR A  14     8468   8562   8647  -1260    -75  -3480       C  
ATOM    107  OH  TYR A  14     -11.881  58.996   7.200  1.00 70.18           O  
ANISOU  107  OH  TYR A  14     8881   9118   8665  -1242   -107  -3697       O  
ATOM    108  N   GLY A  15     -12.829  57.552  17.677  1.00 70.33           N  
ANISOU  108  N   GLY A  15     8426   7852  10442  -1361    311  -1747       N  
ATOM    109  CA  GLY A  15     -12.001  56.473  17.189  1.00 69.48           C  
ANISOU  109  CA  GLY A  15     8374   7462  10560  -1392    366  -1894       C  
ATOM    110  C   GLY A  15     -11.494  56.785  15.792  1.00 67.75           C  
ANISOU  110  C   GLY A  15     8244   7344  10151  -1318    333  -2170       C  
ATOM    111  O   GLY A  15     -10.324  56.588  15.510  1.00 67.50           O  
ANISOU  111  O   GLY A  15     8306   7181  10160  -1230    403  -2214       O  
ATOM    112  N   SER A  16     -14.520  59.746  17.510  1.00 64.73           N  
ANISOU  112  N   SER A  16     7610   7736   9247  -1259    144  -1711       N  
ATOM    113  CA  SER A  16     -13.268  59.872  18.256  1.00 60.88           C  
ANISOU  113  CA  SER A  16     7243   7104   8784  -1159    239  -1519       C  
ATOM    114  C   SER A  16     -12.305  58.751  17.863  1.00 66.60           C  
ANISOU  114  C   SER A  16     8035   7547   9720  -1201    303  -1621       C  
ATOM    115  O   SER A  16     -11.112  58.968  17.689  1.00 64.53           O  
ANISOU  115  O   SER A  16     7889   7216   9412  -1090    345  -1601       O  
ATOM    116  CB  SER A  16     -13.497  59.841  19.762  1.00 58.66           C  
ANISOU  116  CB  SER A  16     6902   6796   8590  -1182    308  -1249       C  
ATOM    117  OG  SER A  16     -14.358  60.904  20.181  1.00 58.58           O  
ANISOU  117  OG  SER A  16     6821   7044   8390  -1124    266  -1175       O  
ATOM    118  N   VAL A  17      -0.741  56.623  -1.910  1.00 92.22           N  
ANISOU  118  N   VAL A  17    12323  12754   9960   -712   1738  -6044       N  
ATOM    119  CA  VAL A  17      -1.060  57.039  -3.284  1.00 97.80           C  
ANISOU  119  CA  VAL A  17    13124  13951  10083   -788   1712  -6239       C  
ATOM    120  C   VAL A  17      -0.238  58.256  -3.688  1.00 96.58           C  
ANISOU  120  C   VAL A  17    13022  14135   9538   -734   1831  -5932       C  
ATOM    121  O   VAL A  17       0.400  58.251  -4.743  1.00100.43           O  
ANISOU  121  O   VAL A  17    13544  14922   9690   -729   2010  -6157       O  
ATOM    122  CB  VAL A  17      -2.585  57.224  -3.546  1.00 97.39           C  
ANISOU  122  CB  VAL A  17    13115  14073   9813   -923   1410  -6256       C  
ATOM    123  CG1 VAL A  17      -2.840  58.010  -4.833  1.00 95.92           C  
ANISOU  123  CG1 VAL A  17    13032  14441   8969   -982   1352  -6288       C  
ATOM    124  CG2 VAL A  17      -3.240  55.856  -3.700  1.00 94.37           C  
ANISOU  124  CG2 VAL A  17    12688  13470   9695  -1012   1350  -6737       C  
ATOM    125  N   ILE A  18       0.124  56.418   0.683  1.00 90.78           N  
ANISOU  125  N   ILE A  18    11988  11815  10687   -550   1784  -5476       N  
ATOM    126  CA  ILE A  18       0.636  55.459  -0.268  1.00 89.62           C  
ANISOU  126  CA  ILE A  18    11827  11678  10545   -533   1962  -5969       C  
ATOM    127  C   ILE A  18       0.389  55.958  -1.674  1.00 92.73           C  
ANISOU  127  C   ILE A  18    12316  12569  10345   -605   1968  -6175       C  
ATOM    128  O   ILE A  18       1.219  55.731  -2.548  1.00100.14           O  
ANISOU  128  O   ILE A  18    13256  13678  11116   -559   2182  -6451       O  
ATOM    129  CB  ILE A  18       0.020  54.059  -0.079  1.00 91.24           C  
ANISOU  129  CB  ILE A  18    11989  11521  11157   -589   1904  -6320       C  
ATOM    130  CG1 ILE A  18       0.347  53.540   1.328  1.00 86.79           C  
ANISOU  130  CG1 ILE A  18    11333  10464  11179   -514   1908  -6077       C  
ATOM    131  CG2 ILE A  18       0.518  53.084  -1.163  1.00 92.84           C  
ANISOU  131  CG2 ILE A  18    12186  11745  11344   -570   2089  -6891       C  
ATOM    132  CD1 ILE A  18       0.220  52.045   1.457  1.00 90.46           C  
ANISOU  132  CD1 ILE A  18    11740  10508  12122   -527   1938  -6434       C  
ATOM    133  N   GLY A  19      -1.443  58.421   1.789  1.00 83.59           N  
ANISOU  133  N   GLY A  19    11148  11110   9501   -629   1388  -4720       N  
ATOM    134  CA  GLY A  19      -0.010  58.556   1.757  1.00 83.29           C  
ANISOU  134  CA  GLY A  19    11087  11051   9509   -521   1614  -4669       C  
ATOM    135  C   GLY A  19       0.613  57.651   0.742  1.00 87.40           C  
ANISOU  135  C   GLY A  19    11593  11616   9996   -506   1808  -5122       C  
ATOM    136  O   GLY A  19       1.522  58.047   0.027  1.00 91.61           O  
ANISOU  136  O   GLY A  19    12141  12383  10283   -460   1983  -5160       O  
ATOM    137  N   ALA A  20      -4.028  57.264   1.493  1.00 82.93           N  
ANISOU  137  N   ALA A  20    11039  10986   9483   -869   1030  -5124       N  
ATOM    138  CA  ALA A  20      -3.666  58.563   0.892  1.00 82.67           C  
ANISOU  138  CA  ALA A  20    11088  11335   8987   -821   1041  -4909       C  
ATOM    139  C   ALA A  20      -2.168  58.728   0.737  1.00 84.53           C  
ANISOU  139  C   ALA A  20    11336  11571   9211   -712   1288  -4874       C  
ATOM    140  O   ALA A  20      -1.677  59.143  -0.306  1.00 93.03           O  
ANISOU  140  O   ALA A  20    12472  12966   9906   -705   1385  -4969       O  
ATOM    141  CB  ALA A  20      -4.188  59.702   1.728  1.00 79.90           C  
ANISOU  141  CB  ALA A  20    10741  11020   8596   -799    885  -4451       C  
ATOM    142  N   SER A  21      -3.959  54.929   3.149  1.00 84.31           N  
ANISOU  142  N   SER A  21    11053  10268  10712   -888   1110  -5307       N  
ATOM    143  CA  SER A  21      -4.144  54.839   1.692  1.00 87.02           C  
ANISOU  143  CA  SER A  21    11454  10928  10679   -948   1113  -5701       C  
ATOM    144  C   SER A  21      -3.663  56.089   0.966  1.00 86.80           C  
ANISOU  144  C   SER A  21    11509  11332  10136   -886   1153  -5548       C  
ATOM    145  O   SER A  21      -2.986  55.978  -0.066  1.00 89.18           O  
ANISOU  145  O   SER A  21    11855  11833  10196   -859   1299  -5818       O  
ATOM    146  CB  SER A  21      -5.602  54.562   1.317  1.00 87.25           C  
ANISOU  146  CB  SER A  21    11469  11063  10616  -1114    893  -5889       C  
ATOM    147  OG  SER A  21      -5.986  53.280   1.754  1.00 85.15           O  
ANISOU  147  OG  SER A  21    11131  10404  10818  -1196    880  -6108       O  
ATOM    148  N   ALA A  22      -3.535  56.148   5.801  1.00 84.81           N  
ANISOU  148  N   ALA A  22    11058  10047  11117   -750   1063  -4375       N  
ATOM    149  CA  ALA A  22      -2.598  55.218   5.170  1.00 85.75           C  
ANISOU  149  CA  ALA A  22    11166  10027  11385   -700   1242  -4709       C  
ATOM    150  C   ALA A  22      -2.760  55.231   3.639  1.00 88.21           C  
ANISOU  150  C   ALA A  22    11540  10663  11311   -751   1264  -5093       C  
ATOM    151  O   ALA A  22      -1.810  55.551   2.908  1.00 88.95           O  
ANISOU  151  O   ALA A  22    11668  10944  11183   -672   1417  -5184       O  
ATOM    152  CB  ALA A  22      -2.820  53.811   5.697  1.00 85.65           C  
ANISOU  152  CB  ALA A  22    11084   9586  11874   -746   1250  -4875       C  
ATOM    153  N   SER A  23      -6.157  57.136   5.909  1.00 79.85           N  
ANISOU  153  N   SER A  23    10417   9768  10154   -941    676  -4199       N  
ATOM    154  CA  SER A  23      -5.050  58.046   5.818  1.00 83.06           C  
ANISOU  154  CA  SER A  23    10889  10305  10365   -812    784  -4009       C  
ATOM    155  C   SER A  23      -3.867  57.306   5.229  1.00 85.48           C  
ANISOU  155  C   SER A  23    11206  10506  10764   -756    983  -4270       C  
ATOM    156  O   SER A  23      -3.281  57.777   4.258  1.00 83.50           O  
ANISOU  156  O   SER A  23    11016  10514  10194   -717   1065  -4372       O  
ATOM    157  CB  SER A  23      -4.772  58.617   7.180  1.00 81.77           C  
ANISOU  157  CB  SER A  23    10700   9991  10376   -740    781  -3602       C  
ATOM    158  OG  SER A  23      -5.040  57.606   8.131  1.00 91.04           O  
ANISOU  158  OG  SER A  23    11797  10827  11968   -784    780  -3590       O  
ATOM    159  N   MET A  24      -7.769  54.771   5.827  1.00 78.19           N  
ANISOU  159  N   MET A  24    10068   9163  10478  -1217    559  -4733       N  
ATOM    160  CA  MET A  24      -8.046  55.947   4.990  1.00 81.68           C  
ANISOU  160  CA  MET A  24    10571  10049  10412  -1191    468  -4698       C  
ATOM    161  C   MET A  24      -6.839  56.880   4.827  1.00 81.93           C  
ANISOU  161  C   MET A  24    10692  10228  10207  -1035    595  -4511       C  
ATOM    162  O   MET A  24      -6.513  57.366   3.733  1.00 89.78           O  
ANISOU  162  O   MET A  24    11761  11522  10827  -1010    620  -4642       O  
ATOM    163  CB  MET A  24      -9.211  56.733   5.574  1.00 81.22           C  
ANISOU  163  CB  MET A  24    10456  10136  10264  -1234    285  -4421       C  
ATOM    164  CG  MET A  24     -10.325  57.017   4.595  1.00 85.37           C  
ANISOU  164  CG  MET A  24    10965  10996  10476  -1332     98  -4604       C  
ATOM    165  SD  MET A  24     -10.834  55.544   3.756  1.00 91.75           S  
ANISOU  165  SD  MET A  24    11725  11697  11438  -1508     64  -5135       S  
ATOM    166  CE  MET A  24     -11.799  56.289   2.465  1.00 96.96           C  
ANISOU  166  CE  MET A  24    12397  12873  11570  -1569   -155  -5282       C  
ATOM    167  N   VAL A  25      -6.824  53.011   7.753  1.00 75.58           N  
ANISOU  167  N   VAL A  25     9643   7996  11077  -1182    745  -4548       N  
ATOM    168  CA  VAL A  25      -6.776  52.613   6.349  1.00 78.65           C  
ANISOU  168  CA  VAL A  25    10072   8518  11291  -1219    767  -5010       C  
ATOM    169  C   VAL A  25      -6.915  53.834   5.435  1.00 79.93           C  
ANISOU  169  C   VAL A  25    10305   9169  10894  -1187    707  -5011       C  
ATOM    170  O   VAL A  25      -6.198  53.953   4.425  1.00 80.99           O  
ANISOU  170  O   VAL A  25    10506   9481  10783  -1126    804  -5239       O  
ATOM    171  CB  VAL A  25      -7.833  51.534   6.048  1.00 79.68           C  
ANISOU  171  CB  VAL A  25    10141   8491  11642  -1409    668  -5329       C  
ATOM    172  CG1 VAL A  25      -7.956  51.259   4.554  1.00 80.06           C  
ANISOU  172  CG1 VAL A  25    10235   8756  11427  -1466    655  -5825       C  
ATOM    173  CG2 VAL A  25      -7.424  50.271   6.769  1.00 80.59           C  
ANISOU  173  CG2 VAL A  25    10207   8090  12323  -1419    762  -5356       C  
ATOM    174  N   ALA A  26      -7.484  54.579   9.963  1.00 74.92           N  
ANISOU  174  N   ALA A  26     9509   8016  10940  -1136    621  -3715       N  
ATOM    175  CA  ALA A  26      -6.092  54.199   9.754  1.00 71.92           C  
ANISOU  175  CA  ALA A  26     9173   7471  10680  -1015    777  -3797       C  
ATOM    176  C   ALA A  26      -5.907  53.803   8.270  1.00 70.77           C  
ANISOU  176  C   ALA A  26     9070   7439  10380  -1034    821  -4249       C  
ATOM    177  O   ALA A  26      -4.995  54.258   7.598  1.00 65.93           O  
ANISOU  177  O   ALA A  26     8515   6986   9547   -927    919  -4325       O  
ATOM    178  CB  ALA A  26      -5.649  53.118  10.737  1.00 68.53           C  
ANISOU  178  CB  ALA A  26     8689   6606  10742  -1009    849  -3690       C  
ATOM    179  N   GLY A  27     -10.265  54.880   9.568  1.00 76.52           N  
ANISOU  179  N   GLY A  27     9582   8571  10921  -1403    317  -3833       N  
ATOM    180  CA  GLY A  27      -9.411  56.032   9.567  1.00 73.50           C  
ANISOU  180  CA  GLY A  27     9291   8376  10259  -1234    364  -3631       C  
ATOM    181  C   GLY A  27      -7.957  55.661   9.340  1.00 75.75           C  
ANISOU  181  C   GLY A  27     9650   8495  10636  -1125    529  -3709       C  
ATOM    182  O   GLY A  27      -7.258  56.344   8.580  1.00 74.64           O  
ANISOU  182  O   GLY A  27     9590   8566  10203  -1030    575  -3760       O  
ATOM    183  N   LEU A  28     -10.888  52.210   9.995  1.00 79.48           N  
ANISOU  183  N   LEU A  28     9818   8278  12100  -1678    353  -4147       N  
ATOM    184  CA  LEU A  28     -11.295  52.880   8.749  1.00 81.36           C  
ANISOU  184  CA  LEU A  28    10090   8905  11917  -1683    253  -4393       C  
ATOM    185  C   LEU A  28     -10.430  54.097   8.506  1.00 78.92           C  
ANISOU  185  C   LEU A  28     9887   8857  11238  -1493    301  -4233       C  
ATOM    186  O   LEU A  28      -9.986  54.334   7.388  1.00 80.63           O  
ANISOU  186  O   LEU A  28    10183   9270  11180  -1445    318  -4468       O  
ATOM    187  CB  LEU A  28     -12.753  53.363   8.798  1.00 81.84           C  
ANISOU  187  CB  LEU A  28    10043   9218  11834  -1803     77  -4332       C  
ATOM    188  CG  LEU A  28     -13.575  53.112   7.524  1.00 85.06           C  
ANISOU  188  CG  LEU A  28    10416   9839  12064  -1934    -63  -4732       C  
ATOM    189  CD1 LEU A  28     -14.561  54.255   7.377  1.00 84.12           C  
ANISOU  189  CD1 LEU A  28    10236  10122  11603  -1928   -230  -4586       C  
ATOM    190  CD2 LEU A  28     -12.741  52.971   6.262  1.00 83.91           C  
ANISOU  190  CD2 LEU A  28    10394   9787  11700  -1873     -9  -5077       C  
ATOM    191  N   LEU A  29      -9.921  51.628  12.646  1.00 71.50           N  
ANISOU  191  N   LEU A  29     8772   6704  11688  -1591    519  -3453       N  
ATOM    192  CA  LEU A  29      -9.333  50.983  11.475  1.00 76.71           C  
ANISOU  192  CA  LEU A  29     9485   7265  12394  -1572    569  -3863       C  
ATOM    193  C   LEU A  29      -9.682  51.674  10.130  1.00 78.29           C  
ANISOU  193  C   LEU A  29     9727   7858  12159  -1574    494  -4163       C  
ATOM    194  O   LEU A  29      -8.854  51.733   9.230  1.00 84.74           O  
ANISOU  194  O   LEU A  29    10623   8746  12825  -1478    565  -4395       O  
ATOM    195  CB  LEU A  29      -9.735  49.504  11.425  1.00 80.04           C  
ANISOU  195  CB  LEU A  29     9844   7299  13266  -1734    571  -4087       C  
ATOM    196  CG  LEU A  29      -8.895  48.384  12.079  1.00 83.37           C  
ANISOU  196  CG  LEU A  29    10265   7227  14183  -1696    681  -4008       C  
ATOM    197  CD1 LEU A  29      -7.906  48.836  13.133  1.00 80.60           C  
ANISOU  197  CD1 LEU A  29     9944   6817  13861  -1522    758  -3590       C  
ATOM    198  CD2 LEU A  29      -9.816  47.323  12.663  1.00 81.93           C  
ANISOU  198  CD2 LEU A  29     9987   6724  14417  -1908    635  -3965       C  
ATOM    199  N   PHE A  30     -11.893  53.345  13.999  1.00 74.27           N  
ANISOU  199  N   PHE A  30     8986   7559  11672  -1673    352  -2939       N  
ATOM    200  CA  PHE A  30     -10.435  53.512  14.126  1.00 71.09           C  
ANISOU  200  CA  PHE A  30     8697   7047  11264  -1497    457  -2848       C  
ATOM    201  C   PHE A  30      -9.670  52.899  12.933  1.00 71.27           C  
ANISOU  201  C   PHE A  30     8787   6960  11331  -1463    508  -3215       C  
ATOM    202  O   PHE A  30      -8.873  53.580  12.288  1.00 68.66           O  
ANISOU  202  O   PHE A  30     8543   6802  10742  -1328    544  -3285       O  
ATOM    203  CB  PHE A  30      -9.948  52.873  15.427  1.00 69.28           C  
ANISOU  203  CB  PHE A  30     8448   6504  11369  -1492    537  -2547       C  
ATOM    204  CG  PHE A  30      -8.463  52.721  15.503  1.00 70.80           C  
ANISOU  204  CG  PHE A  30     8724   6519  11658  -1332    635  -2498       C  
ATOM    205  CD1 PHE A  30      -7.659  53.806  15.806  1.00 72.37           C  
ANISOU  205  CD1 PHE A  30     8987   6903  11607  -1169    657  -2301       C  
ATOM    206  CD2 PHE A  30      -7.855  51.515  15.234  1.00 74.91           C  
ANISOU  206  CD2 PHE A  30     9245   6684  12531  -1341    703  -2667       C  
ATOM    207  CE1 PHE A  30      -6.283  53.680  15.868  1.00 70.36           C  
ANISOU  207  CE1 PHE A  30     8782   6501  11449  -1026    743  -2257       C  
ATOM    208  CE2 PHE A  30      -6.481  51.383  15.280  1.00 73.42           C  
ANISOU  208  CE2 PHE A  30     9107   6342  12445  -1179    794  -2630       C  
ATOM    209  CZ  PHE A  30      -5.696  52.469  15.600  1.00 72.26           C  
ANISOU  209  CZ  PHE A  30     9009   6405  12040  -1025    813  -2420       C  
ATOM    210  N   ALA A  31     -14.263  51.971  13.343  1.00 79.28           N  
ANISOU  210  N   ALA A  31     9380   8145  12597  -2085    186  -3327       N  
ATOM    211  CA  ALA A  31     -14.060  53.330  12.895  1.00 75.44           C  
ANISOU  211  CA  ALA A  31     8965   8031  11668  -1920    139  -3291       C  
ATOM    212  C   ALA A  31     -12.565  53.670  12.892  1.00 75.27           C  
ANISOU  212  C   ALA A  31     9097   7939  11563  -1721    249  -3212       C  
ATOM    213  O   ALA A  31     -12.038  54.210  11.900  1.00 76.31           O  
ANISOU  213  O   ALA A  31     9321   8256  11417  -1620    239  -3387       O  
ATOM    214  CB  ALA A  31     -14.840  54.279  13.775  1.00 70.59           C  
ANISOU  214  CB  ALA A  31     8265   7650  10906  -1899     95  -2983       C  
ATOM    215  N   VAL A  32     -13.771  49.624  14.718  1.00 81.52           N  
ANISOU  215  N   VAL A  32     9616   7596  13759  -2270    348  -3149       N  
ATOM    216  CA  VAL A  32     -13.901  49.565  13.278  1.00 83.38           C  
ANISOU  216  CA  VAL A  32     9879   7954  13845  -2292    276  -3601       C  
ATOM    217  C   VAL A  32     -13.629  50.972  12.753  1.00 81.26           C  
ANISOU  217  C   VAL A  32     9688   8113  13073  -2112    236  -3592       C  
ATOM    218  O   VAL A  32     -12.785  51.164  11.877  1.00 87.42           O  
ANISOU  218  O   VAL A  32    10580   8952  13683  -1990    266  -3802       O  
ATOM    219  CB  VAL A  32     -15.266  48.953  12.847  1.00 83.45           C  
ANISOU  219  CB  VAL A  32     9743   7974  13989  -2548    161  -3838       C  
ATOM    220  CG1 VAL A  32     -15.551  49.176  11.370  1.00 80.96           C  
ANISOU  220  CG1 VAL A  32     9449   7916  13395  -2565     49  -4279       C  
ATOM    221  CG2 VAL A  32     -15.264  47.472  13.158  1.00 80.73           C  
ANISOU  221  CG2 VAL A  32     9357   7136  14178  -2715    213  -3904       C  
TER     222      VAL A  32                                                      
ENDMDL                                                                          
MODEL        2                                                                  
ATOM      1  N   GLY A   1       8.058  59.006  12.665  1.00 80.38           N  
ATOM      2  CA  GLY A   1       8.986  57.987  12.141  1.00 83.37           C  
ATOM      3  C   GLY A   1       9.621  58.375  10.813  1.00 84.33           C  
ATOM      4  O   GLY A   1      10.857  58.412  10.667  1.00 72.14           O  
ATOM      5  N   VAL A   2       6.399  61.381  12.666  1.00 84.79           N  
ATOM      6  CA  VAL A   2       7.486  61.188  13.644  1.00 83.84           C  
ATOM      7  C   VAL A   2       8.506  60.185  13.086  1.00 82.59           C  
ATOM      8  O   VAL A   2       9.690  60.499  13.007  1.00 80.71           O  
ATOM      9  CB  VAL A   2       7.005  60.800  15.065  1.00 81.10           C  
ATOM     10  CG1 VAL A   2       8.135  60.959  16.081  1.00 81.48           C  
ATOM     11  CG2 VAL A   2       5.843  61.679  15.489  1.00 80.73           C  
ATOM     12  N   ARG A   3       4.787  61.058  10.259  1.00 98.19           N  
ATOM     13  CA  ARG A   3       5.462  62.335  10.588  1.00 94.68           C  
ATOM     14  C   ARG A   3       6.542  62.260  11.678  1.00 90.57           C  
ATOM     15  O   ARG A   3       7.503  63.029  11.617  1.00 95.24           O  
ATOM     16  CB  ARG A   3       4.433  63.416  10.920  1.00 97.47           C  
ATOM     17  CG  ARG A   3       4.962  64.669  11.631  1.00100.52           C  
ATOM     18  CD  ARG A   3       4.522  65.935  10.933  1.00 99.14           C  
ATOM     19  NE  ARG A   3       5.346  66.163   9.750  1.00 98.23           N  
ATOM     20  CZ  ARG A   3       5.097  67.074   8.819  1.00101.11           C  
ATOM     21  NH1 ARG A   3       4.031  67.878   8.907  1.00100.32           N  
ATOM     22  NH2 ARG A   3       5.932  67.184   7.790  1.00104.85           N  
ATOM     23  N   GLY A   4       3.940  58.236  10.061  1.00 98.28           N  
ATOM     24  CA  GLY A   4       4.240  59.092   8.895  1.00100.67           C  
ATOM     25  C   GLY A   4       5.056  60.365   9.148  1.00103.54           C  
ATOM     26  O   GLY A   4       5.935  60.702   8.349  1.00112.56           O  
ATOM     27  N   PHE A   5       3.404  57.550  12.915  1.00 85.39           N  
ATOM     28  CA  PHE A   5       4.372  56.881  12.035  1.00 89.96           C  
ATOM     29  C   PHE A   5       4.868  57.794  10.915  1.00 93.64           C  
ATOM     30  O   PHE A   5       6.066  58.073  10.832  1.00 94.85           O  
ATOM     31  CB  PHE A   5       3.792  55.605  11.413  1.00 95.55           C  
ATOM     32  CG  PHE A   5       3.792  54.406  12.331  1.00 97.35           C  
ATOM     33  CD1 PHE A   5       4.960  53.977  12.955  1.00 97.44           C  
ATOM     34  CD2 PHE A   5       2.628  53.675  12.533  1.00 93.99           C  
ATOM     35  CE1 PHE A   5       4.956  52.867  13.779  1.00 94.40           C  
ATOM     36  CE2 PHE A   5       2.624  52.559  13.343  1.00 91.55           C  
ATOM     37  CZ  PHE A   5       3.789  52.156  13.969  1.00 93.60           C  
ATOM     38  N   LEU A   6       1.500  59.606  13.693  1.00 86.86           N  
ATOM     39  CA  LEU A   6       2.680  59.273  14.521  1.00 87.89           C  
ATOM     40  C   LEU A   6       3.783  58.508  13.760  1.00 87.20           C  
ATOM     41  O   LEU A   6       4.965  58.775  13.962  1.00 84.49           O  
ATOM     42  CB  LEU A   6       2.259  58.462  15.751  1.00 90.68           C  
ATOM     43  CG  LEU A   6       3.317  58.039  16.780  1.00 92.60           C  
ATOM     44  CD1 LEU A   6       3.927  59.236  17.490  1.00 95.59           C  
ATOM     45  CD2 LEU A   6       2.691  57.095  17.790  1.00 91.88           C  
ATOM     46  N   ALA A   7      -0.503  59.130  11.667  1.00 88.27           N  
ATOM     47  CA  ALA A   7       0.320  60.338  11.651  1.00 86.96           C  
ATOM     48  C   ALA A   7       1.589  60.245  12.522  1.00 86.52           C  
ATOM     49  O   ALA A   7       2.635  60.742  12.115  1.00 86.56           O  
ATOM     50  CB  ALA A   7      -0.513  61.547  12.041  1.00 90.39           C  
ATOM     51  N   MET A   8      -1.759  56.997  13.061  1.00 88.85           N  
ATOM     52  CA  MET A   8      -0.936  56.773  11.879  1.00 91.83           C  
ATOM     53  C   MET A   8       0.028  57.910  11.669  1.00 91.48           C  
ATOM     54  O   MET A   8       1.228  57.681  11.516  1.00 89.90           O  
ATOM     55  CB  MET A   8      -1.800  56.611  10.635  1.00 97.68           C  
ATOM     56  CG  MET A   8      -2.401  55.234  10.488  1.00100.85           C  
ATOM     57  SD  MET A   8      -1.129  53.963  10.579  1.00105.73           S  
ATOM     58  CE  MET A   8      -1.443  53.337  12.233  1.00118.11           C  
ATOM     59  N   SER A   9      -3.104  58.390  15.089  1.00 76.97           N  
ATOM     60  CA  SER A   9      -2.224  57.257  15.420  1.00 80.80           C  
ATOM     61  C   SER A   9      -1.249  56.959  14.290  1.00 84.65           C  
ATOM     62  O   SER A   9      -0.069  56.694  14.533  1.00 80.30           O  
ATOM     63  CB  SER A   9      -3.021  55.986  15.739  1.00 83.11           C  
ATOM     64  OG  SER A   9      -3.159  55.802  17.137  1.00 88.07           O  
ATOM     65  N   SER A  10      -4.659  60.035  13.347  1.00 78.44           N  
ATOM     66  CA  SER A  10      -3.647  60.607  14.243  1.00 78.36           C  
ATOM     67  C   SER A  10      -2.624  59.560  14.670  1.00 76.40           C  
ATOM     68  O   SER A  10      -1.428  59.821  14.646  1.00 77.82           O  
ATOM     69  CB  SER A  10      -4.310  61.193  15.496  1.00 80.79           C  
ATOM     70  OG  SER A  10      -3.416  62.024  16.222  1.00 82.99           O  
ATOM     71  N   GLY A  11      -6.392  58.056  12.300  1.00 80.31           N  
ATOM     72  CA  GLY A  11      -5.555  58.884  11.424  1.00 78.65           C  
ATOM     73  C   GLY A  11      -4.385  59.539  12.140  1.00 80.24           C  
ATOM     74  O   GLY A  11      -3.268  59.607  11.599  1.00 81.32           O  
ATOM     75  N   GLY A  12      -7.617  57.158  14.794  1.00 79.52           N  
ATOM     76  CA  GLY A  12      -6.969  56.260  13.822  1.00 83.16           C  
ATOM     77  C   GLY A  12      -5.967  56.929  12.882  1.00 83.52           C  
ATOM     78  O   GLY A  12      -4.845  56.424  12.690  1.00 80.28           O  
ATOM     79  N   LEU A  13      -8.829  59.390  16.057  1.00 80.77           N  
ATOM     80  CA  LEU A  13      -7.674  58.730  16.671  1.00 78.56           C  
ATOM     81  C   LEU A  13      -6.912  57.819  15.707  1.00 81.26           C  
ATOM     82  O   LEU A  13      -5.691  57.712  15.806  1.00 87.75           O  
ATOM     83  CB  LEU A  13      -8.096  57.930  17.909  1.00 77.34           C  
ATOM     84  CG  LEU A  13      -7.046  57.063  18.633  1.00 74.05           C  
ATOM     85  CD1 LEU A  13      -6.008  57.885  19.384  1.00 71.04           C  
ATOM     86  CD2 LEU A  13      -7.751  56.116  19.582  1.00 74.50           C  
ATOM     87  N   GLY A  14     -10.938  59.775  14.113  1.00 85.16           N  
ATOM     88  CA  GLY A  14      -9.992  60.824  14.486  1.00 85.86           C  
ATOM     89  C   GLY A  14      -8.704  60.226  15.027  1.00 87.77           C  
ATOM     90  O   GLY A  14      -7.622  60.494  14.496  1.00 94.35           O  
ATOM     91  N   TYR A  15     -12.188  57.250  14.520  1.00 74.32           N  
ATOM     92  CA  TYR A  15     -11.660  57.640  13.196  1.00 78.10           C  
ATOM     93  C   TYR A  15     -10.614  58.752  13.320  1.00 80.86           C  
ATOM     94  O   TYR A  15      -9.544  58.673  12.717  1.00 80.60           O  
ATOM     95  CB  TYR A  15     -12.823  58.101  12.304  1.00 78.18           C  
ATOM     96  CG  TYR A  15     -12.496  58.397  10.858  1.00 77.60           C  
ATOM     97  CD1 TYR A  15     -12.397  57.371   9.919  1.00 79.83           C  
ATOM     98  CD2 TYR A  15     -12.339  59.712  10.415  1.00 79.66           C  
ATOM     99  CE1 TYR A  15     -12.119  57.640   8.583  1.00 81.89           C  
ATOM    100  CE2 TYR A  15     -12.061  59.995   9.083  1.00 84.26           C  
ATOM    101  CZ  TYR A  15     -11.951  58.957   8.169  1.00 85.66           C  
ATOM    102  OH  TYR A  15     -11.676  59.245   6.850  1.00 86.51           O  
ATOM    103  N   GLY A  16     -13.313  57.509  17.131  1.00 67.46           N  
ATOM    104  CA  GLY A  16     -12.142  56.652  16.874  1.00 67.58           C  
ATOM    105  C   GLY A  16     -11.418  56.905  15.551  1.00 69.32           C  
ATOM    106  O   GLY A  16     -10.187  56.777  15.461  1.00 63.60           O  
ATOM    107  N   THR A  17     -15.707  59.079  16.716  1.00 78.72           N  
ATOM    108  CA  THR A  17     -14.532  59.651  17.408  1.00 76.19           C  
ATOM    109  C   THR A  17     -13.211  58.825  17.319  1.00 71.43           C  
ATOM    110  O   THR A  17     -12.123  59.397  17.439  1.00 70.66           O  
ATOM    111  CB  THR A  17     -14.877  60.042  18.883  1.00 74.40           C  
ATOM    112  OG1 THR A  17     -14.586  61.427  19.095  1.00 81.86           O  
ATOM    113  CG2 THR A  17     -14.121  59.225  19.909  1.00 73.57           C  
ATOM    114  N   GLY A  18      -0.307  59.327  -2.737  1.00 67.33           N  
ATOM    115  CA  GLY A  18       0.496  60.540  -2.776  1.00 67.32           C  
ATOM    116  C   GLY A  18       1.975  60.281  -2.974  1.00 68.67           C  
ATOM    117  O   GLY A  18       2.629  60.972  -3.754  1.00 70.56           O  
ATOM    118  N   LEU A  19      -0.992  56.651  -2.261  1.00 71.55           N  
ATOM    119  CA  LEU A  19      -1.326  57.287  -3.542  1.00 71.39           C  
ATOM    120  C   LEU A  19      -0.524  58.558  -3.803  1.00 70.39           C  
ATOM    121  O   LEU A  19      -0.119  58.819  -4.940  1.00 71.13           O  
ATOM    122  CB  LEU A  19      -2.812  57.621  -3.590  1.00 72.00           C  
ATOM    123  CG  LEU A  19      -3.699  56.393  -3.771  1.00 72.75           C  
ATOM    124  CD1 LEU A  19      -5.133  56.652  -3.316  1.00 73.94           C  
ATOM    125  CD2 LEU A  19      -3.660  55.922  -5.219  1.00 71.01           C  
ATOM    126  N   PHE A  20      -0.104  56.693   0.468  1.00 86.19           N  
ATOM    127  CA  PHE A  20       0.440  55.756  -0.534  1.00 78.57           C  
ATOM    128  C   PHE A  20       0.252  56.304  -1.934  1.00 76.29           C  
ATOM    129  O   PHE A  20       1.210  56.384  -2.700  1.00 80.88           O  
ATOM    130  CB  PHE A  20      -0.235  54.383  -0.447  1.00 78.63           C  
ATOM    131  CG  PHE A  20       0.284  53.379  -1.447  1.00 77.91           C  
ATOM    132  CD1 PHE A  20       1.370  52.565  -1.137  1.00 76.59           C  
ATOM    133  CD2 PHE A  20      -0.320  53.240  -2.695  1.00 76.30           C  
ATOM    134  CE1 PHE A  20       1.860  51.642  -2.044  1.00 75.90           C  
ATOM    135  CE2 PHE A  20       0.161  52.317  -3.610  1.00 77.40           C  
ATOM    136  CZ  PHE A  20       1.254  51.516  -3.283  1.00 78.33           C  
ATOM    137  N   GLY A  21      -1.842  58.649   1.587  1.00 84.05           N  
ATOM    138  CA  GLY A  21      -0.383  58.829   1.618  1.00 87.23           C  
ATOM    139  C   GLY A  21       0.391  57.921   0.664  1.00 90.58           C  
ATOM    140  O   GLY A  21       1.438  58.325   0.134  1.00 94.29           O  
ATOM    141  N   SER A  22      -4.232  57.284   1.260  1.00 86.55           N  
ATOM    142  CA  SER A  22      -4.063  58.580   0.586  1.00 85.23           C  
ATOM    143  C   SER A  22      -2.606  59.036   0.562  1.00 84.66           C  
ATOM    144  O   SER A  22      -2.197  59.718  -0.373  1.00 89.02           O  
ATOM    145  CB  SER A  22      -4.941  59.652   1.219  1.00 85.10           C  
ATOM    146  OG  SER A  22      -6.308  59.312   1.092  1.00 83.17           O  
ATOM    147  N   VAL A  23      -3.693  54.944   2.969  1.00 75.30           N  
ATOM    148  CA  VAL A  23      -3.905  54.875   1.529  1.00 79.24           C  
ATOM    149  C   VAL A  23      -3.611  56.194   0.816  1.00 84.91           C  
ATOM    150  O   VAL A  23      -2.818  56.220  -0.125  1.00 88.99           O  
ATOM    151  CB  VAL A  23      -5.328  54.375   1.208  1.00 81.24           C  
ATOM    152  CG1 VAL A  23      -5.752  54.708  -0.218  1.00 82.02           C  
ATOM    153  CG2 VAL A  23      -5.387  52.876   1.430  1.00 83.26           C  
ATOM    154  N   LEU A  24      -3.241  55.970   5.751  1.00 78.32           N  
ATOM    155  CA  LEU A  24      -2.388  55.028   4.998  1.00 80.32           C  
ATOM    156  C   LEU A  24      -2.499  55.183   3.488  1.00 78.29           C  
ATOM    157  O   LEU A  24      -1.523  55.538   2.825  1.00 81.62           O  
ATOM    158  CB  LEU A  24      -2.696  53.566   5.386  1.00 80.86           C  
ATOM    159  CG  LEU A  24      -1.828  52.827   6.426  1.00 82.69           C  
ATOM    160  CD1 LEU A  24      -0.757  51.981   5.753  1.00 84.24           C  
ATOM    161  CD2 LEU A  24      -1.186  53.753   7.439  1.00 80.64           C  
ATOM    162  N   SER A  25      -5.653  57.311   6.146  1.00 75.64           N  
ATOM    163  CA  SER A  25      -4.367  58.045   6.269  1.00 75.50           C  
ATOM    164  C   SER A  25      -3.355  57.264   5.442  1.00 77.16           C  
ATOM    165  O   SER A  25      -2.743  57.812   4.521  1.00 79.43           O  
ATOM    166  CB  SER A  25      -3.895  58.185   7.747  1.00 74.94           C  
ATOM    167  OG  SER A  25      -3.204  59.411   8.041  1.00 68.33           O  
ATOM    168  N   THR A  26      -7.439  55.056   5.824  1.00 65.86           N  
ATOM    169  CA  THR A  26      -7.523  56.225   4.945  1.00 73.62           C  
ATOM    170  C   THR A  26      -6.256  57.116   4.967  1.00 74.88           C  
ATOM    171  O   THR A  26      -5.843  57.617   3.914  1.00 75.74           O  
ATOM    172  CB  THR A  26      -8.766  57.065   5.264  1.00 77.48           C  
ATOM    173  OG1 THR A  26      -9.860  56.185   5.545  1.00 82.41           O  
ATOM    174  CG2 THR A  26      -9.130  57.987   4.087  1.00 77.38           C  
ATOM    175  N   VAL A  27      -6.580  53.403   7.935  1.00 65.93           N  
ATOM    176  CA  VAL A  27      -6.464  52.939   6.562  1.00 64.22           C  
ATOM    177  C   VAL A  27      -6.558  54.105   5.577  1.00 65.79           C  
ATOM    178  O   VAL A  27      -5.809  54.154   4.602  1.00 72.15           O  
ATOM    179  CB  VAL A  27      -7.464  51.793   6.287  1.00 61.00           C  
ATOM    180  CG1 VAL A  27      -7.616  51.488   4.798  1.00 59.09           C  
ATOM    181  CG2 VAL A  27      -6.993  50.539   7.020  1.00 61.26           C  
ATOM    182  N   GLY A  28      -7.146  54.938  10.214  1.00 75.61           N  
ATOM    183  CA  GLY A  28      -5.786  54.427  10.004  1.00 75.21           C  
ATOM    184  C   GLY A  28      -5.518  53.781   8.645  1.00 73.11           C  
ATOM    185  O   GLY A  28      -4.360  53.601   8.258  1.00 75.35           O  
ATOM    186  N   GLY A  29      -9.995  54.900   9.697  1.00 80.60           N  
ATOM    187  CA  GLY A  29      -9.176  56.114   9.693  1.00 79.36           C  
ATOM    188  C   GLY A  29      -7.740  55.780   9.372  1.00 77.71           C  
ATOM    189  O   GLY A  29      -7.207  56.223   8.354  1.00 79.80           O  
ATOM    190  N   LEU A  30     -10.680  52.136   9.974  1.00 68.15           N  
ATOM    191  CA  LEU A  30     -11.040  52.912   8.797  1.00 70.17           C  
ATOM    192  C   LEU A  30     -10.146  54.136   8.619  1.00 77.31           C  
ATOM    193  O   LEU A  30      -9.592  54.365   7.534  1.00 81.98           O  
ATOM    194  CB  LEU A  30     -12.486  53.366   8.885  1.00 67.59           C  
ATOM    195  CG  LEU A  30     -13.554  52.289   8.912  1.00 64.44           C  
ATOM    196  CD1 LEU A  30     -14.915  52.953   9.008  1.00 65.48           C  
ATOM    197  CD2 LEU A  30     -13.460  51.424   7.675  1.00 62.92           C  
ATOM    198  N   LEU A  31      -9.791  51.504  12.571  1.00 84.78           N  
ATOM    199  CA  LEU A  31      -9.188  50.876  11.395  1.00 82.80           C  
ATOM    200  C   LEU A  31      -9.444  51.691  10.145  1.00 76.45           C  
ATOM    201  O   LEU A  31      -8.530  51.913   9.355  1.00 76.25           O  
ATOM    202  CB  LEU A  31      -9.717  49.466  11.185  1.00 81.53           C  
ATOM    203  CG  LEU A  31      -9.517  48.516  12.357  1.00 79.91           C  
ATOM    204  CD1 LEU A  31     -10.083  47.176  11.921  1.00 81.70           C  
ATOM    205  CD2 LEU A  31      -8.060  48.408  12.814  1.00 77.06           C  
ATOM    206  N   PHE A  32     -11.485  53.249  14.058  1.00 88.10           N  
ATOM    207  CA  PHE A  32     -10.033  53.239  14.263  1.00 88.48           C  
ATOM    208  C   PHE A  32      -9.260  52.593  13.117  1.00 86.25           C  
ATOM    209  O   PHE A  32      -8.188  53.067  12.755  1.00 90.90           O  
ATOM    210  CB  PHE A  32      -9.664  52.518  15.556  1.00 92.14           C  
ATOM    211  CG  PHE A  32      -8.183  52.437  15.786  1.00 95.28           C  
ATOM    212  CD1 PHE A  32      -7.479  53.551  16.232  1.00 98.94           C  
ATOM    213  CD2 PHE A  32      -7.484  51.268  15.518  1.00 97.81           C  
ATOM    214  CE1 PHE A  32      -6.108  53.496  16.433  1.00102.50           C  
ATOM    215  CE2 PHE A  32      -6.116  51.206  15.717  1.00105.01           C  
ATOM    216  CZ  PHE A  32      -5.424  52.324  16.167  1.00105.34           C  
ATOM    217  N   SER A  33     -13.968  52.120  13.342  1.00 91.90           N  
ATOM    218  CA  SER A  33     -13.558  53.413  12.773  1.00 91.90           C  
ATOM    219  C   SER A  33     -12.043  53.633  12.909  1.00 89.11           C  
ATOM    220  O   SER A  33     -11.398  54.100  11.970  1.00 83.38           O  
ATOM    221  CB  SER A  33     -14.366  54.575  13.382  1.00 88.96           C  
ATOM    222  OG  SER A  33     -13.965  54.881  14.704  1.00 85.11           O  
ATOM    223  N   VAL A  34     -13.896  49.748  14.975  1.00 87.47           N  
ATOM    224  CA  VAL A  34     -14.087  49.721  13.520  1.00 86.78           C  
ATOM    225  C   VAL A  34     -13.504  50.959  12.879  1.00 89.36           C  
ATOM    226  O   VAL A  34     -12.679  50.860  11.969  1.00 90.48           O  
ATOM    227  CB  VAL A  34     -15.569  49.668  13.136  1.00 86.66           C  
ATOM    228  CG1 VAL A  34     -15.770  49.936  11.644  1.00 82.52           C  
ATOM    229  CG2 VAL A  34     -16.135  48.315  13.518  1.00 93.14           C  
ATOM    230  N   THR A  35     -13.550  50.695  17.591  1.00 89.35           N  
ATOM    231  CA  THR A  35     -12.513  49.840  17.010  1.00 90.77           C  
ATOM    232  C   THR A  35     -12.710  49.495  15.525  1.00 89.57           C  
ATOM    233  O   THR A  35     -11.787  49.001  14.888  1.00 90.71           O  
ATOM    234  CB  THR A  35     -12.382  48.551  17.827  1.00 92.76           C  
ATOM    235  OG1 THR A  35     -12.373  48.891  19.220  1.00 95.65           O  
ATOM    236  CG2 THR A  35     -11.094  47.795  17.484  1.00 93.50           C  
ATOM    237  N   SER A  36     -16.091  51.946  17.865  1.00 84.15           N  
ATOM    238  CA  SER A  36     -14.839  52.725  17.945  1.00 83.59           C  
ATOM    239  C   SER A  36     -13.671  51.977  17.276  1.00 85.46           C  
ATOM    240  O   SER A  36     -12.921  52.545  16.480  1.00 82.35           O  
ATOM    241  CB  SER A  36     -14.481  53.094  19.392  1.00 78.68           C  
ATOM    242  OG  SER A  36     -15.054  54.342  19.759  1.00 73.19           O  
TER     243      SER A  36                                                      
ENDMDL                                                                          
END                                                                             """

ENSEMBLE_POLYALA = """MODEL        1                                                                  
ATOM      1  N   ALA A   1       6.808  61.263  13.432  1.00 72.61           N  
ANISOU    1  N   ALA A   1     9095   7963  10530    149   1518  -1669       N  
ATOM      2  CA  ALA A   1       7.802  60.901  14.433  1.00 75.34           C  
ANISOU    2  CA  ALA A   1     9317   8124  11185    233   1511  -1496       C  
ATOM      3  C   ALA A   1       8.710  59.751  13.957  1.00 83.54           C  
ANISOU    3  C   ALA A   1    10223   8997  12519    336   1668  -1701       C  
ATOM      4  O   ALA A   1       9.847  59.594  14.417  1.00 89.53           O  
ANISOU    4  O   ALA A   1    10845   9666  13503    427   1706  -1596       O  
ATOM      5  CB  ALA A   1       7.145  60.691  15.789  1.00 75.55           C  
ANISOU    5  CB  ALA A   1     9364   7997  11342    213   1344  -1273       C  
ATOM      6  N   ALA A   2       5.296  61.322  10.847  1.00 71.18           N  
ANISOU    6  N   ALA A   2     9096   8147   9801     20   1606  -2181       N  
ATOM      7  CA  ALA A   2       6.152  62.390  11.390  1.00 71.60           C  
ANISOU    7  CA  ALA A   2     9122   8276   9806     47   1598  -1902       C  
ATOM      8  C   ALA A   2       7.216  61.908  12.350  1.00 76.14           C  
ANISOU    8  C   ALA A   2     9562   8635  10730    135   1627  -1774       C  
ATOM      9  O   ALA A   2       8.425  62.127  12.148  1.00 78.15           O  
ANISOU    9  O   ALA A   2     9721   8935  11036    188   1743  -1750       O  
ATOM     10  CB  ALA A   2       5.322  63.389  12.190  1.00 81.86           C  
ANISOU   10  CB  ALA A   2    10508   9630  10963     -5   1414  -1647       C  
ATOM     11  N   ALA A   3       4.070  58.691  10.381  1.00 78.09           N  
ANISOU   11  N   ALA A   3     9947   8649  11074    -14   1625  -2713       N  
ATOM     12  CA  ALA A   3       4.632  59.476   9.274  1.00 77.89           C  
ANISOU   12  CA  ALA A   3     9952   8925  10717     -7   1737  -2805       C  
ATOM     13  C   ALA A   3       5.589  60.662   9.714  1.00 79.72           C  
ANISOU   13  C   ALA A   3    10158   9278  10852     34   1755  -2495       C  
ATOM     14  O   ALA A   3       6.580  60.956   9.019  1.00 79.47           O  
ANISOU   14  O   ALA A   3    10082   9389  10722     72   1910  -2556       O  
ATOM     15  CB  ALA A   3       3.483  59.989   8.383  1.00 78.46           C  
ANISOU   15  CB  ALA A   3    10144   9263  10404   -112   1654  -2920       C  
ATOM     16  N   ALA A   4       3.150  58.127  13.054  1.00 73.86           N  
ANISOU   16  N   ALA A   4     9370   7685  11008    -44   1372  -2230       N  
ATOM     17  CA  ALA A   4       4.164  57.312  12.379  1.00 79.41           C  
ANISOU   17  CA  ALA A   4     9998   8267  11906     42   1533  -2465       C  
ATOM     18  C   ALA A   4       4.849  58.044  11.235  1.00 79.34           C  
ANISOU   18  C   ALA A   4    10005   8539  11602     70   1658  -2607       C  
ATOM     19  O   ALA A   4       6.064  58.036  11.132  1.00 82.59           O  
ANISOU   19  O   ALA A   4    10329   8935  12116    166   1782  -2605       O  
ATOM     20  CB  ALA A   4       3.529  56.054  11.792  1.00 90.18           C  
ANISOU   20  CB  ALA A   4    11362   9448  13451     -3   1563  -2793       C  
ATOM     21  N   ALA A   5       1.211  60.152  13.195  1.00 76.46           N  
ANISOU   21  N   ALA A   5     9861   8451  10738   -188   1133  -1999       N  
ATOM     22  CA  ALA A   5       2.272  60.073  14.171  1.00 72.86           C  
ANISOU   22  CA  ALA A   5     9339   7838  10506   -107   1163  -1788       C  
ATOM     23  C   ALA A   5       3.430  59.310  13.580  1.00 72.45           C  
ANISOU   23  C   ALA A   5     9214   7668  10645    -28   1316  -1965       C  
ATOM     24  O   ALA A   5       4.556  59.785  13.621  1.00 73.23           O  
ANISOU   24  O   ALA A   5     9266   7815  10740     44   1382  -1871       O  
ATOM     25  CB  ALA A   5       1.830  59.378  15.430  1.00 75.99           C  
ANISOU   25  CB  ALA A   5     9696   8009  11167   -128   1083  -1618       C  
ATOM     26  N   ALA A   6      -0.524  59.587  11.161  1.00 66.38           N  
ANISOU   26  N   ALA A   6     8660   7395   9166   -337   1099  -2523       N  
ATOM     27  CA  ALA A   6       0.233  60.855  11.098  1.00 63.60           C  
ANISOU   27  CA  ALA A   6     8346   7236   8581   -274   1122  -2331       C  
ATOM     28  C   ALA A   6       1.411  60.765  12.032  1.00 67.49           C  
ANISOU   28  C   ALA A   6     8775   7548   9320   -188   1188  -2138       C  
ATOM     29  O   ALA A   6       2.505  61.188  11.690  1.00 69.06           O  
ANISOU   29  O   ALA A   6     8959   7820   9461   -127   1289  -2114       O  
ATOM     30  CB  ALA A   6      -0.614  62.070  11.372  1.00 56.93           C  
ANISOU   30  CB  ALA A   6     7562   6584   7483   -304    980  -2128       C  
ATOM     31  N   ALA A   7      -1.601  57.282  11.899  1.00 61.88           N  
ANISOU   31  N   ALA A   7     7992   6336   9182   -453   1058  -2717       N  
ATOM     32  CA  ALA A   7      -0.766  57.188  10.752  1.00 67.84           C  
ANISOU   32  CA  ALA A   7     8764   7177   9832   -397   1183  -2972       C  
ATOM     33  C   ALA A   7       0.013  58.494  10.624  1.00 70.11           C  
ANISOU   33  C   ALA A   7     9092   7710   9835   -319   1216  -2784       C  
ATOM     34  O   ALA A   7       1.136  58.503  10.080  1.00 69.35           O  
ANISOU   34  O   ALA A   7     8979   7645   9722   -241   1355  -2874       O  
ATOM     35  CB  ALA A   7      -1.441  56.665   9.450  1.00 73.01           C  
ANISOU   35  CB  ALA A   7     9448   7939  10350   -483   1185  -3372       C  
ATOM     36  N   ALA A   8      -2.856  58.098  14.278  1.00 60.88           N  
ANISOU   36  N   ALA A   8     7836   6188   9105   -513    856  -2123       N  
ATOM     37  CA  ALA A   8      -1.791  57.132  14.351  1.00 59.72           C  
ANISOU   37  CA  ALA A   8     7651   5777   9263   -452    963  -2180       C  
ATOM     38  C   ALA A   8      -0.982  57.170  13.066  1.00 60.98           C  
ANISOU   38  C   ALA A   8     7835   6031   9301   -391   1071  -2449       C  
ATOM     39  O   ALA A   8       0.223  57.129  13.131  1.00 62.29           O  
ANISOU   39  O   ALA A   8     7974   6123   9569   -287   1166  -2409       O  
ATOM     40  CB  ALA A   8      -2.352  55.726  14.631  1.00 62.50           C  
ANISOU   40  CB  ALA A   8     7947   5828   9972   -539    960  -2276       C  
ATOM     41  N   ALA A   9      -4.846  59.990  13.395  1.00 58.40           N  
ANISOU   41  N   ALA A   9     7589   6436   8161   -595    652  -2149       N  
ATOM     42  CA  ALA A   9      -3.712  60.360  14.217  1.00 60.43           C  
ANISOU   42  CA  ALA A   9     7859   6604   8495   -498    713  -1918       C  
ATOM     43  C   ALA A   9      -2.570  59.381  14.182  1.00 58.45           C  
ANISOU   43  C   ALA A   9     7580   6117   8509   -453    834  -1997       C  
ATOM     44  O   ALA A   9      -1.430  59.784  14.089  1.00 59.36           O  
ANISOU   44  O   ALA A   9     7710   6250   8590   -363    905  -1939       O  
ATOM     45  CB  ALA A   9      -4.092  60.532  15.675  1.00 58.59           C  
ANISOU   45  CB  ALA A   9     7590   6300   8370   -505    657  -1640       C  
ATOM     46  N   ALA A  10      -6.665  58.234  12.008  1.00 64.10           N  
ANISOU   46  N   ALA A  10     8232   7113   9010   -830    573  -2697       N  
ATOM     47  CA  ALA A  10      -5.894  59.264  11.295  1.00 66.28           C  
ANISOU   47  CA  ALA A  10     8594   7613   8975   -726    604  -2683       C  
ATOM     48  C   ALA A  10      -4.683  59.763  12.098  1.00 63.40           C  
ANISOU   48  C   ALA A  10     8255   7161   8671   -608    688  -2424       C  
ATOM     49  O   ALA A  10      -3.609  59.922  11.539  1.00 68.98           O  
ANISOU   49  O   ALA A  10     9001   7896   9309   -536    786  -2482       O  
ATOM     50  CB  ALA A  10      -6.774  60.431  10.945  1.00 62.41           C  
ANISOU   50  CB  ALA A  10     8127   7431   8153   -729    472  -2610       C  
ATOM     51  N   ALA A  11      -7.523  56.981  14.399  1.00 63.45           N  
ANISOU   51  N   ALA A  11     8001   6600   9506   -957    560  -2345       N  
ATOM     52  CA  ALA A  11      -6.926  56.234  13.306  1.00 64.48           C  
ANISOU   52  CA  ALA A  11     8161   6628   9707   -955    625  -2661       C  
ATOM     53  C   ALA A  11      -6.072  57.139  12.464  1.00 63.66           C  
ANISOU   53  C   ALA A  11     8144   6743   9299   -835    664  -2722       C  
ATOM     54  O   ALA A  11      -4.886  56.872  12.245  1.00 60.49           O  
ANISOU   54  O   ALA A  11     7771   6229   8982   -747    777  -2775       O  
ATOM     55  N   ALA A  12      -8.666  58.912  15.975  1.00 61.98           N  
ANISOU   55  N   ALA A  12     7779   6763   9007   -914    436  -1866       N  
ATOM     56  CA  ALA A  12      -7.453  58.270  16.446  1.00 60.60           C  
ANISOU   56  CA  ALA A  12     7634   6335   9053   -862    535  -1781       C  
ATOM     57  C   ALA A  12      -6.762  57.549  15.313  1.00 64.32           C  
ANISOU   57  C   ALA A  12     8139   6692   9606   -853    596  -2064       C  
ATOM     58  O   ALA A  12      -5.548  57.496  15.275  1.00 65.65           O  
ANISOU   58  O   ALA A  12     8348   6762   9831   -754    676  -2042       O  
ATOM     59  CB  ALA A  12      -7.713  57.288  17.567  1.00 60.13           C  
ANISOU   59  CB  ALA A  12     7504   6036   9306   -947    562  -1617       C  
ATOM     60  N   ALA A  13     -10.929  59.546  14.435  1.00 68.55           N  
ANISOU   60  N   ALA A  13     8519   8026   9501  -1036    226  -2184       N  
ATOM     61  CA  ALA A  13      -9.869  60.510  14.575  1.00 63.92           C  
ANISOU   61  CA  ALA A  13     8040   7496   8748   -885    267  -2029       C  
ATOM     62  C   ALA A  13      -8.582  59.839  15.040  1.00 63.57           C  
ANISOU   62  C   ALA A  13     8042   7192   8919   -846    392  -1965       C  
ATOM     63  O   ALA A  13      -7.516  60.131  14.527  1.00 72.56           O  
ANISOU   63  O   ALA A  13     9262   8339   9969   -758    447  -2000       O  
ATOM     64  CB  ALA A  13     -10.261  61.598  15.526  1.00 62.67           C  
ANISOU   64  CB  ALA A  13     7865   7469   8476   -815    225  -1776       C  
ATOM     65  N   ALA A  14     -12.378  57.236  14.902  1.00 66.70           N  
ANISOU   65  N   ALA A  14     8073   7456   9813  -1356    227  -2355       N  
ATOM     66  CA  ALA A  14     -11.949  57.558  13.541  1.00 67.37           C  
ANISOU   66  CA  ALA A  14     8246   7680   9670  -1296    192  -2604       C  
ATOM     67  C   ALA A  14     -10.808  58.532  13.580  1.00 65.65           C  
ANISOU   67  C   ALA A  14     8152   7540   9251  -1118    245  -2455       C  
ATOM     68  O   ALA A  14      -9.829  58.360  12.875  1.00 64.93           O  
ANISOU   68  O   ALA A  14     8148   7399   9123  -1062    310  -2590       O  
ATOM     69  CB  ALA A  14     -13.052  58.190  12.723  1.00 66.46           C  
ANISOU   69  CB  ALA A  14     8075   7875   9300  -1332     45  -2738       C  
ATOM     70  N   ALA A  15     -12.829  57.552  17.677  1.00 70.33           N  
ANISOU   70  N   ALA A  15     8426   7852  10442  -1361    311  -1747       N  
ATOM     71  CA  ALA A  15     -12.001  56.473  17.189  1.00 69.48           C  
ANISOU   71  CA  ALA A  15     8374   7462  10560  -1392    366  -1894       C  
ATOM     72  C   ALA A  15     -11.494  56.785  15.792  1.00 67.75           C  
ANISOU   72  C   ALA A  15     8244   7344  10151  -1318    333  -2170       C  
ATOM     73  O   ALA A  15     -10.324  56.588  15.510  1.00 67.50           O  
ANISOU   73  O   ALA A  15     8306   7181  10160  -1230    403  -2214       O  
ATOM     74  N   ALA A  16     -14.520  59.746  17.510  1.00 64.73           N  
ANISOU   74  N   ALA A  16     7610   7736   9247  -1259    144  -1711       N  
ATOM     75  CA  ALA A  16     -13.268  59.872  18.256  1.00 60.88           C  
ANISOU   75  CA  ALA A  16     7243   7104   8784  -1159    239  -1519       C  
ATOM     76  C   ALA A  16     -12.305  58.751  17.863  1.00 66.60           C  
ANISOU   76  C   ALA A  16     8035   7547   9720  -1201    303  -1621       C  
ATOM     77  O   ALA A  16     -11.112  58.968  17.689  1.00 64.53           O  
ANISOU   77  O   ALA A  16     7889   7216   9412  -1090    345  -1601       O  
ATOM     78  CB  ALA A  16     -13.497  59.841  19.762  1.00 58.66           C  
ANISOU   78  CB  ALA A  16     6902   6796   8590  -1182    308  -1249       C  
ATOM     79  N   ALA A  17      -0.741  56.623  -1.910  1.00 92.22           N  
ANISOU   79  N   ALA A  17    12323  12754   9960   -712   1738  -6044       N  
ATOM     80  CA  ALA A  17      -1.060  57.039  -3.284  1.00 97.80           C  
ANISOU   80  CA  ALA A  17    13124  13951  10083   -788   1712  -6239       C  
ATOM     81  C   ALA A  17      -0.238  58.256  -3.688  1.00 96.58           C  
ANISOU   81  C   ALA A  17    13022  14135   9538   -734   1831  -5932       C  
ATOM     82  O   ALA A  17       0.400  58.251  -4.743  1.00100.43           O  
ANISOU   82  O   ALA A  17    13544  14922   9690   -729   2010  -6157       O  
ATOM     83  CB  ALA A  17      -2.585  57.224  -3.546  1.00 97.39           C  
ANISOU   83  CB  ALA A  17    13115  14073   9813   -923   1410  -6256       C  
ATOM     84  N   ALA A  18       0.124  56.418   0.683  1.00 90.78           N  
ANISOU   84  N   ALA A  18    11988  11815  10687   -550   1784  -5476       N  
ATOM     85  CA  ALA A  18       0.636  55.459  -0.268  1.00 89.62           C  
ANISOU   85  CA  ALA A  18    11827  11678  10545   -533   1962  -5969       C  
ATOM     86  C   ALA A  18       0.389  55.958  -1.674  1.00 92.73           C  
ANISOU   86  C   ALA A  18    12316  12569  10345   -605   1968  -6175       C  
ATOM     87  O   ALA A  18       1.219  55.731  -2.548  1.00100.14           O  
ANISOU   87  O   ALA A  18    13256  13678  11116   -559   2182  -6451       O  
ATOM     88  CB  ALA A  18       0.020  54.059  -0.079  1.00 91.24           C  
ANISOU   88  CB  ALA A  18    11989  11521  11157   -589   1904  -6320       C  
ATOM     89  N   ALA A  19      -1.443  58.421   1.789  1.00 83.59           N  
ANISOU   89  N   ALA A  19    11148  11110   9501   -629   1388  -4720       N  
ATOM     90  CA  ALA A  19      -0.010  58.556   1.757  1.00 83.29           C  
ANISOU   90  CA  ALA A  19    11087  11051   9509   -521   1614  -4669       C  
ATOM     91  C   ALA A  19       0.613  57.651   0.742  1.00 87.40           C  
ANISOU   91  C   ALA A  19    11593  11616   9996   -506   1808  -5122       C  
ATOM     92  O   ALA A  19       1.522  58.047   0.027  1.00 91.61           O  
ANISOU   92  O   ALA A  19    12141  12383  10283   -460   1983  -5160       O  
ATOM     93  N   ALA A  20      -4.028  57.264   1.493  1.00 82.93           N  
ANISOU   93  N   ALA A  20    11039  10986   9483   -869   1030  -5124       N  
ATOM     94  CA  ALA A  20      -3.666  58.563   0.892  1.00 82.67           C  
ANISOU   94  CA  ALA A  20    11088  11335   8987   -821   1041  -4909       C  
ATOM     95  C   ALA A  20      -2.168  58.728   0.737  1.00 84.53           C  
ANISOU   95  C   ALA A  20    11336  11571   9211   -712   1288  -4874       C  
ATOM     96  O   ALA A  20      -1.677  59.143  -0.306  1.00 93.03           O  
ANISOU   96  O   ALA A  20    12472  12966   9906   -705   1385  -4969       O  
ATOM     97  CB  ALA A  20      -4.188  59.702   1.728  1.00 79.90           C  
ANISOU   97  CB  ALA A  20    10741  11020   8596   -799    885  -4451       C  
ATOM     98  N   ALA A  21      -3.959  54.929   3.149  1.00 84.31           N  
ANISOU   98  N   ALA A  21    11053  10268  10712   -888   1110  -5307       N  
ATOM     99  CA  ALA A  21      -4.144  54.839   1.692  1.00 87.02           C  
ANISOU   99  CA  ALA A  21    11454  10928  10679   -948   1113  -5701       C  
ATOM    100  C   ALA A  21      -3.663  56.089   0.966  1.00 86.80           C  
ANISOU  100  C   ALA A  21    11509  11332  10136   -886   1153  -5548       C  
ATOM    101  O   ALA A  21      -2.986  55.978  -0.066  1.00 89.18           O  
ANISOU  101  O   ALA A  21    11855  11833  10196   -859   1299  -5818       O  
ATOM    102  CB  ALA A  21      -5.602  54.562   1.317  1.00 87.25           C  
ANISOU  102  CB  ALA A  21    11469  11063  10616  -1114    893  -5889       C  
ATOM    103  N   ALA A  22      -3.535  56.148   5.801  1.00 84.81           N  
ANISOU  103  N   ALA A  22    11058  10047  11117   -750   1063  -4375       N  
ATOM    104  CA  ALA A  22      -2.598  55.218   5.170  1.00 85.75           C  
ANISOU  104  CA  ALA A  22    11166  10027  11385   -700   1242  -4709       C  
ATOM    105  C   ALA A  22      -2.760  55.231   3.639  1.00 88.21           C  
ANISOU  105  C   ALA A  22    11540  10663  11311   -751   1264  -5093       C  
ATOM    106  O   ALA A  22      -1.810  55.551   2.908  1.00 88.95           O  
ANISOU  106  O   ALA A  22    11668  10944  11183   -672   1417  -5184       O  
ATOM    107  CB  ALA A  22      -2.820  53.811   5.697  1.00 85.65           C  
ANISOU  107  CB  ALA A  22    11084   9586  11874   -746   1250  -4875       C  
ATOM    108  N   ALA A  23      -6.157  57.136   5.909  1.00 79.85           N  
ANISOU  108  N   ALA A  23    10417   9768  10154   -941    676  -4199       N  
ATOM    109  CA  ALA A  23      -5.050  58.046   5.818  1.00 83.06           C  
ANISOU  109  CA  ALA A  23    10889  10305  10365   -812    784  -4009       C  
ATOM    110  C   ALA A  23      -3.867  57.306   5.229  1.00 85.48           C  
ANISOU  110  C   ALA A  23    11206  10506  10764   -756    983  -4270       C  
ATOM    111  O   ALA A  23      -3.281  57.777   4.258  1.00 83.50           O  
ANISOU  111  O   ALA A  23    11016  10514  10194   -717   1065  -4372       O  
ATOM    112  CB  ALA A  23      -4.772  58.617   7.180  1.00 81.77           C  
ANISOU  112  CB  ALA A  23    10700   9991  10376   -740    781  -3602       C  
ATOM    113  N   ALA A  24      -7.769  54.771   5.827  1.00 78.19           N  
ANISOU  113  N   ALA A  24    10068   9163  10478  -1217    559  -4733       N  
ATOM    114  CA  ALA A  24      -8.046  55.947   4.990  1.00 81.68           C  
ANISOU  114  CA  ALA A  24    10571  10049  10412  -1191    468  -4698       C  
ATOM    115  C   ALA A  24      -6.839  56.880   4.827  1.00 81.93           C  
ANISOU  115  C   ALA A  24    10692  10228  10207  -1035    595  -4511       C  
ATOM    116  O   ALA A  24      -6.513  57.366   3.733  1.00 89.78           O  
ANISOU  116  O   ALA A  24    11761  11522  10827  -1010    620  -4642       O  
ATOM    117  CB  ALA A  24      -9.211  56.733   5.574  1.00 81.22           C  
ANISOU  117  CB  ALA A  24    10456  10136  10264  -1234    285  -4421       C  
ATOM    118  N   ALA A  25      -6.824  53.011   7.753  1.00 75.58           N  
ANISOU  118  N   ALA A  25     9643   7996  11077  -1182    745  -4548       N  
ATOM    119  CA  ALA A  25      -6.776  52.613   6.349  1.00 78.65           C  
ANISOU  119  CA  ALA A  25    10072   8518  11291  -1219    767  -5010       C  
ATOM    120  C   ALA A  25      -6.915  53.834   5.435  1.00 79.93           C  
ANISOU  120  C   ALA A  25    10305   9169  10894  -1187    707  -5011       C  
ATOM    121  O   ALA A  25      -6.198  53.953   4.425  1.00 80.99           O  
ANISOU  121  O   ALA A  25    10506   9481  10783  -1126    804  -5239       O  
ATOM    122  CB  ALA A  25      -7.833  51.534   6.048  1.00 79.68           C  
ANISOU  122  CB  ALA A  25    10141   8491  11642  -1409    668  -5329       C  
ATOM    123  N   ALA A  26      -7.484  54.579   9.963  1.00 74.92           N  
ANISOU  123  N   ALA A  26     9509   8016  10940  -1136    621  -3715       N  
ATOM    124  CA  ALA A  26      -6.092  54.199   9.754  1.00 71.92           C  
ANISOU  124  CA  ALA A  26     9173   7471  10680  -1015    777  -3797       C  
ATOM    125  C   ALA A  26      -5.907  53.803   8.270  1.00 70.77           C  
ANISOU  125  C   ALA A  26     9070   7439  10380  -1034    821  -4249       C  
ATOM    126  O   ALA A  26      -4.995  54.258   7.598  1.00 65.93           O  
ANISOU  126  O   ALA A  26     8515   6986   9547   -927    919  -4325       O  
ATOM    127  CB  ALA A  26      -5.649  53.118  10.737  1.00 68.53           C  
ANISOU  127  CB  ALA A  26     8689   6606  10742  -1009    849  -3690       C  
ATOM    128  N   ALA A  27     -10.265  54.880   9.568  1.00 76.52           N  
ANISOU  128  N   ALA A  27     9582   8571  10921  -1403    317  -3833       N  
ATOM    129  CA  ALA A  27      -9.411  56.032   9.567  1.00 73.50           C  
ANISOU  129  CA  ALA A  27     9291   8376  10259  -1234    364  -3631       C  
ATOM    130  C   ALA A  27      -7.957  55.661   9.340  1.00 75.75           C  
ANISOU  130  C   ALA A  27     9650   8495  10636  -1125    529  -3709       C  
ATOM    131  O   ALA A  27      -7.258  56.344   8.580  1.00 74.64           O  
ANISOU  131  O   ALA A  27     9590   8566  10203  -1030    575  -3760       O  
ATOM    132  N   ALA A  28     -10.888  52.210   9.995  1.00 79.48           N  
ANISOU  132  N   ALA A  28     9818   8278  12100  -1678    353  -4147       N  
ATOM    133  CA  ALA A  28     -11.295  52.880   8.749  1.00 81.36           C  
ANISOU  133  CA  ALA A  28    10090   8905  11917  -1683    253  -4393       C  
ATOM    134  C   ALA A  28     -10.430  54.097   8.506  1.00 78.92           C  
ANISOU  134  C   ALA A  28     9887   8857  11238  -1493    301  -4233       C  
ATOM    135  O   ALA A  28      -9.986  54.334   7.388  1.00 80.63           O  
ANISOU  135  O   ALA A  28    10183   9270  11180  -1445    318  -4468       O  
ATOM    136  CB  ALA A  28     -12.753  53.363   8.798  1.00 81.84           C  
ANISOU  136  CB  ALA A  28    10043   9218  11834  -1803     77  -4332       C  
ATOM    137  N   ALA A  29      -9.921  51.628  12.646  1.00 71.50           N  
ANISOU  137  N   ALA A  29     8772   6704  11688  -1591    519  -3453       N  
ATOM    138  CA  ALA A  29      -9.333  50.983  11.475  1.00 76.71           C  
ANISOU  138  CA  ALA A  29     9485   7265  12394  -1572    569  -3863       C  
ATOM    139  C   ALA A  29      -9.682  51.674  10.130  1.00 78.29           C  
ANISOU  139  C   ALA A  29     9727   7858  12159  -1574    494  -4163       C  
ATOM    140  O   ALA A  29      -8.854  51.733   9.230  1.00 84.74           O  
ANISOU  140  O   ALA A  29    10623   8746  12825  -1478    565  -4395       O  
ATOM    141  CB  ALA A  29      -9.735  49.504  11.425  1.00 80.04           C  
ANISOU  141  CB  ALA A  29     9844   7299  13266  -1734    571  -4087       C  
ATOM    142  N   ALA A  30     -11.893  53.345  13.999  1.00 74.27           N  
ANISOU  142  N   ALA A  30     8986   7559  11672  -1673    352  -2939       N  
ATOM    143  CA  ALA A  30     -10.435  53.512  14.126  1.00 71.09           C  
ANISOU  143  CA  ALA A  30     8697   7047  11264  -1497    457  -2848       C  
ATOM    144  C   ALA A  30      -9.670  52.899  12.933  1.00 71.27           C  
ANISOU  144  C   ALA A  30     8787   6960  11331  -1463    508  -3215       C  
ATOM    145  O   ALA A  30      -8.873  53.580  12.288  1.00 68.66           O  
ANISOU  145  O   ALA A  30     8543   6802  10742  -1328    544  -3285       O  
ATOM    146  CB  ALA A  30      -9.948  52.873  15.427  1.00 69.28           C  
ANISOU  146  CB  ALA A  30     8448   6504  11369  -1492    537  -2547       C  
ATOM    147  N   ALA A  31     -14.263  51.971  13.343  1.00 79.28           N  
ANISOU  147  N   ALA A  31     9380   8145  12597  -2085    186  -3327       N  
ATOM    148  CA  ALA A  31     -14.060  53.330  12.895  1.00 75.44           C  
ANISOU  148  CA  ALA A  31     8965   8031  11668  -1920    139  -3291       C  
ATOM    149  C   ALA A  31     -12.565  53.670  12.892  1.00 75.27           C  
ANISOU  149  C   ALA A  31     9097   7939  11563  -1721    249  -3212       C  
ATOM    150  O   ALA A  31     -12.038  54.210  11.900  1.00 76.31           O  
ANISOU  150  O   ALA A  31     9321   8256  11417  -1620    239  -3387       O  
ATOM    151  CB  ALA A  31     -14.840  54.279  13.775  1.00 70.59           C  
ANISOU  151  CB  ALA A  31     8265   7650  10906  -1899     95  -2983       C  
ATOM    152  N   ALA A  32     -13.771  49.624  14.718  1.00 81.52           N  
ANISOU  152  N   ALA A  32     9616   7596  13759  -2270    348  -3149       N  
ATOM    153  CA  ALA A  32     -13.901  49.565  13.278  1.00 83.38           C  
ANISOU  153  CA  ALA A  32     9879   7954  13845  -2292    276  -3601       C  
ATOM    154  C   ALA A  32     -13.629  50.972  12.753  1.00 81.26           C  
ANISOU  154  C   ALA A  32     9688   8113  13073  -2112    236  -3592       C  
ATOM    155  O   ALA A  32     -12.785  51.164  11.877  1.00 87.42           O  
ANISOU  155  O   ALA A  32    10580   8952  13683  -1990    266  -3802       O  
ATOM    156  CB  ALA A  32     -15.266  48.953  12.847  1.00 83.45           C  
ANISOU  156  CB  ALA A  32     9743   7974  13989  -2548    161  -3838       C  
TER     157      ALA A  32                                                      
ENDMDL                                                                          
MODEL        2                                                                  
ATOM      1  N   ALA A   1       8.058  59.006  12.665  1.00 80.38           N  
ATOM      2  CA  ALA A   1       8.986  57.987  12.141  1.00 83.37           C  
ATOM      3  C   ALA A   1       9.621  58.375  10.813  1.00 84.33           C  
ATOM      4  O   ALA A   1      10.857  58.412  10.667  1.00 72.14           O  
ATOM      5  N   ALA A   2       6.399  61.381  12.666  1.00 84.79           N  
ATOM      6  CA  ALA A   2       7.486  61.188  13.644  1.00 83.84           C  
ATOM      7  C   ALA A   2       8.506  60.185  13.086  1.00 82.59           C  
ATOM      8  O   ALA A   2       9.690  60.499  13.007  1.00 80.71           O  
ATOM      9  CB  ALA A   2       7.005  60.800  15.065  1.00 81.10           C  
ATOM     10  N   ALA A   3       4.787  61.058  10.259  1.00 98.19           N  
ATOM     11  CA  ALA A   3       5.462  62.335  10.588  1.00 94.68           C  
ATOM     12  C   ALA A   3       6.542  62.260  11.678  1.00 90.57           C  
ATOM     13  O   ALA A   3       7.503  63.029  11.617  1.00 95.24           O  
ATOM     14  CB  ALA A   3       4.433  63.416  10.920  1.00 97.47           C  
ATOM     15  N   ALA A   4       3.940  58.236  10.061  1.00 98.28           N  
ATOM     16  CA  ALA A   4       4.240  59.092   8.895  1.00100.67           C  
ATOM     17  C   ALA A   4       5.056  60.365   9.148  1.00103.54           C  
ATOM     18  O   ALA A   4       5.935  60.702   8.349  1.00112.56           O  
ATOM     19  N   ALA A   5       3.404  57.550  12.915  1.00 85.39           N  
ATOM     20  CA  ALA A   5       4.372  56.881  12.035  1.00 89.96           C  
ATOM     21  C   ALA A   5       4.868  57.794  10.915  1.00 93.64           C  
ATOM     22  O   ALA A   5       6.066  58.073  10.832  1.00 94.85           O  
ATOM     23  CB  ALA A   5       3.792  55.605  11.413  1.00 95.55           C  
ATOM     24  N   ALA A   6       1.500  59.606  13.693  1.00 86.86           N  
ATOM     25  CA  ALA A   6       2.680  59.273  14.521  1.00 87.89           C  
ATOM     26  C   ALA A   6       3.783  58.508  13.760  1.00 87.20           C  
ATOM     27  O   ALA A   6       4.965  58.775  13.962  1.00 84.49           O  
ATOM     28  CB  ALA A   6       2.259  58.462  15.751  1.00 90.68           C  
ATOM     29  N   ALA A   7      -0.503  59.130  11.667  1.00 88.27           N  
ATOM     30  CA  ALA A   7       0.320  60.338  11.651  1.00 86.96           C  
ATOM     31  C   ALA A   7       1.589  60.245  12.522  1.00 86.52           C  
ATOM     32  O   ALA A   7       2.635  60.742  12.115  1.00 86.56           O  
ATOM     33  CB  ALA A   7      -0.513  61.547  12.041  1.00 90.39           C  
ATOM     34  N   ALA A   8      -1.759  56.997  13.061  1.00 88.85           N  
ATOM     35  CA  ALA A   8      -0.936  56.773  11.879  1.00 91.83           C  
ATOM     36  C   ALA A   8       0.028  57.910  11.669  1.00 91.48           C  
ATOM     37  O   ALA A   8       1.228  57.681  11.516  1.00 89.90           O  
ATOM     38  CB  ALA A   8      -1.800  56.611  10.635  1.00 97.68           C  
ATOM     39  N   ALA A   9      -3.104  58.390  15.089  1.00 76.97           N  
ATOM     40  CA  ALA A   9      -2.224  57.257  15.420  1.00 80.80           C  
ATOM     41  C   ALA A   9      -1.249  56.959  14.290  1.00 84.65           C  
ATOM     42  O   ALA A   9      -0.069  56.694  14.533  1.00 80.30           O  
ATOM     43  CB  ALA A   9      -3.021  55.986  15.739  1.00 83.11           C  
ATOM     44  N   ALA A  10      -4.659  60.035  13.347  1.00 78.44           N  
ATOM     45  CA  ALA A  10      -3.647  60.607  14.243  1.00 78.36           C  
ATOM     46  C   ALA A  10      -2.624  59.560  14.670  1.00 76.40           C  
ATOM     47  O   ALA A  10      -1.428  59.821  14.646  1.00 77.82           O  
ATOM     48  CB  ALA A  10      -4.310  61.193  15.496  1.00 80.79           C  
ATOM     49  N   ALA A  11      -6.392  58.056  12.300  1.00 80.31           N  
ATOM     50  CA  ALA A  11      -5.555  58.884  11.424  1.00 78.65           C  
ATOM     51  C   ALA A  11      -4.385  59.539  12.140  1.00 80.24           C  
ATOM     52  O   ALA A  11      -3.268  59.607  11.599  1.00 81.32           O  
ATOM     53  N   ALA A  12      -7.617  57.158  14.794  1.00 79.52           N  
ATOM     54  CA  ALA A  12      -6.969  56.260  13.822  1.00 83.16           C  
ATOM     55  C   ALA A  12      -5.967  56.929  12.882  1.00 83.52           C  
ATOM     56  O   ALA A  12      -4.845  56.424  12.690  1.00 80.28           O  
ATOM     57  N   ALA A  13      -8.829  59.390  16.057  1.00 80.77           N  
ATOM     58  CA  ALA A  13      -7.674  58.730  16.671  1.00 78.56           C  
ATOM     59  C   ALA A  13      -6.912  57.819  15.707  1.00 81.26           C  
ATOM     60  O   ALA A  13      -5.691  57.712  15.806  1.00 87.75           O  
ATOM     61  CB  ALA A  13      -8.096  57.930  17.909  1.00 77.34           C  
ATOM     62  N   ALA A  14     -10.938  59.775  14.113  1.00 85.16           N  
ATOM     63  CA  ALA A  14      -9.992  60.824  14.486  1.00 85.86           C  
ATOM     64  C   ALA A  14      -8.704  60.226  15.027  1.00 87.77           C  
ATOM     65  O   ALA A  14      -7.622  60.494  14.496  1.00 94.35           O  
ATOM     66  N   ALA A  15     -12.188  57.250  14.520  1.00 74.32           N  
ATOM     67  CA  ALA A  15     -11.660  57.640  13.196  1.00 78.10           C  
ATOM     68  C   ALA A  15     -10.614  58.752  13.320  1.00 80.86           C  
ATOM     69  O   ALA A  15      -9.544  58.673  12.717  1.00 80.60           O  
ATOM     70  CB  ALA A  15     -12.823  58.101  12.304  1.00 78.18           C  
ATOM     71  N   ALA A  16     -13.313  57.509  17.131  1.00 67.46           N  
ATOM     72  CA  ALA A  16     -12.142  56.652  16.874  1.00 67.58           C  
ATOM     73  C   ALA A  16     -11.418  56.905  15.551  1.00 69.32           C  
ATOM     74  O   ALA A  16     -10.187  56.777  15.461  1.00 63.60           O  
ATOM     75  N   ALA A  17     -15.707  59.079  16.716  1.00 78.72           N  
ATOM     76  CA  ALA A  17     -14.532  59.651  17.408  1.00 76.19           C  
ATOM     77  C   ALA A  17     -13.211  58.825  17.319  1.00 71.43           C  
ATOM     78  O   ALA A  17     -12.123  59.397  17.439  1.00 70.66           O  
ATOM     79  CB  ALA A  17     -14.877  60.042  18.883  1.00 74.40           C  
ATOM     80  N   ALA A  18      -0.307  59.327  -2.737  1.00 67.33           N  
ATOM     81  CA  ALA A  18       0.496  60.540  -2.776  1.00 67.32           C  
ATOM     82  C   ALA A  18       1.975  60.281  -2.974  1.00 68.67           C  
ATOM     83  O   ALA A  18       2.629  60.972  -3.754  1.00 70.56           O  
ATOM     84  N   ALA A  19      -0.992  56.651  -2.261  1.00 71.55           N  
ATOM     85  CA  ALA A  19      -1.326  57.287  -3.542  1.00 71.39           C  
ATOM     86  C   ALA A  19      -0.524  58.558  -3.803  1.00 70.39           C  
ATOM     87  O   ALA A  19      -0.119  58.819  -4.940  1.00 71.13           O  
ATOM     88  CB  ALA A  19      -2.812  57.621  -3.590  1.00 72.00           C  
ATOM     89  N   ALA A  20      -0.104  56.693   0.468  1.00 86.19           N  
ATOM     90  CA  ALA A  20       0.440  55.756  -0.534  1.00 78.57           C  
ATOM     91  C   ALA A  20       0.252  56.304  -1.934  1.00 76.29           C  
ATOM     92  O   ALA A  20       1.210  56.384  -2.700  1.00 80.88           O  
ATOM     93  CB  ALA A  20      -0.235  54.383  -0.447  1.00 78.63           C  
ATOM     94  N   ALA A  21      -1.842  58.649   1.587  1.00 84.05           N  
ATOM     95  CA  ALA A  21      -0.383  58.829   1.618  1.00 87.23           C  
ATOM     96  C   ALA A  21       0.391  57.921   0.664  1.00 90.58           C  
ATOM     97  O   ALA A  21       1.438  58.325   0.134  1.00 94.29           O  
ATOM     98  N   ALA A  22      -4.232  57.284   1.260  1.00 86.55           N  
ATOM     99  CA  ALA A  22      -4.063  58.580   0.586  1.00 85.23           C  
ATOM    100  C   ALA A  22      -2.606  59.036   0.562  1.00 84.66           C  
ATOM    101  O   ALA A  22      -2.197  59.718  -0.373  1.00 89.02           O  
ATOM    102  CB  ALA A  22      -4.941  59.652   1.219  1.00 85.10           C  
ATOM    103  N   ALA A  23      -3.693  54.944   2.969  1.00 75.30           N  
ATOM    104  CA  ALA A  23      -3.905  54.875   1.529  1.00 79.24           C  
ATOM    105  C   ALA A  23      -3.611  56.194   0.816  1.00 84.91           C  
ATOM    106  O   ALA A  23      -2.818  56.220  -0.125  1.00 88.99           O  
ATOM    107  CB  ALA A  23      -5.328  54.375   1.208  1.00 81.24           C  
ATOM    108  N   ALA A  24      -3.241  55.970   5.751  1.00 78.32           N  
ATOM    109  CA  ALA A  24      -2.388  55.028   4.998  1.00 80.32           C  
ATOM    110  C   ALA A  24      -2.499  55.183   3.488  1.00 78.29           C  
ATOM    111  O   ALA A  24      -1.523  55.538   2.825  1.00 81.62           O  
ATOM    112  CB  ALA A  24      -2.696  53.566   5.386  1.00 80.86           C  
ATOM    113  N   ALA A  25      -5.653  57.311   6.146  1.00 75.64           N  
ATOM    114  CA  ALA A  25      -4.367  58.045   6.269  1.00 75.50           C  
ATOM    115  C   ALA A  25      -3.355  57.264   5.442  1.00 77.16           C  
ATOM    116  O   ALA A  25      -2.743  57.812   4.521  1.00 79.43           O  
ATOM    117  CB  ALA A  25      -3.895  58.185   7.747  1.00 74.94           C  
ATOM    118  N   ALA A  26      -7.439  55.056   5.824  1.00 65.86           N  
ATOM    119  CA  ALA A  26      -7.523  56.225   4.945  1.00 73.62           C  
ATOM    120  C   ALA A  26      -6.256  57.116   4.967  1.00 74.88           C  
ATOM    121  O   ALA A  26      -5.843  57.617   3.914  1.00 75.74           O  
ATOM    122  CB  ALA A  26      -8.766  57.065   5.264  1.00 77.48           C  
ATOM    123  N   ALA A  27      -6.580  53.403   7.935  1.00 65.93           N  
ATOM    124  CA  ALA A  27      -6.464  52.939   6.562  1.00 64.22           C  
ATOM    125  C   ALA A  27      -6.558  54.105   5.577  1.00 65.79           C  
ATOM    126  O   ALA A  27      -5.809  54.154   4.602  1.00 72.15           O  
ATOM    127  CB  ALA A  27      -7.464  51.793   6.287  1.00 61.00           C  
ATOM    128  N   ALA A  28      -7.146  54.938  10.214  1.00 75.61           N  
ATOM    129  CA  ALA A  28      -5.786  54.427  10.004  1.00 75.21           C  
ATOM    130  C   ALA A  28      -5.518  53.781   8.645  1.00 73.11           C  
ATOM    131  O   ALA A  28      -4.360  53.601   8.258  1.00 75.35           O  
ATOM    132  N   ALA A  29      -9.995  54.900   9.697  1.00 80.60           N  
ATOM    133  CA  ALA A  29      -9.176  56.114   9.693  1.00 79.36           C  
ATOM    134  C   ALA A  29      -7.740  55.780   9.372  1.00 77.71           C  
ATOM    135  O   ALA A  29      -7.207  56.223   8.354  1.00 79.80           O  
ATOM    136  N   ALA A  30     -10.680  52.136   9.974  1.00 68.15           N  
ATOM    137  CA  ALA A  30     -11.040  52.912   8.797  1.00 70.17           C  
ATOM    138  C   ALA A  30     -10.146  54.136   8.619  1.00 77.31           C  
ATOM    139  O   ALA A  30      -9.592  54.365   7.534  1.00 81.98           O  
ATOM    140  CB  ALA A  30     -12.486  53.366   8.885  1.00 67.59           C  
ATOM    141  N   ALA A  31      -9.791  51.504  12.571  1.00 84.78           N  
ATOM    142  CA  ALA A  31      -9.188  50.876  11.395  1.00 82.80           C  
ATOM    143  C   ALA A  31      -9.444  51.691  10.145  1.00 76.45           C  
ATOM    144  O   ALA A  31      -8.530  51.913   9.355  1.00 76.25           O  
ATOM    145  CB  ALA A  31      -9.717  49.466  11.185  1.00 81.53           C  
ATOM    146  N   ALA A  32     -11.485  53.249  14.058  1.00 88.10           N  
ATOM    147  CA  ALA A  32     -10.033  53.239  14.263  1.00 88.48           C  
ATOM    148  C   ALA A  32      -9.260  52.593  13.117  1.00 86.25           C  
ATOM    149  O   ALA A  32      -8.188  53.067  12.755  1.00 90.90           O  
ATOM    150  CB  ALA A  32      -9.664  52.518  15.556  1.00 92.14           C  
ATOM    151  N   ALA A  33     -13.968  52.120  13.342  1.00 91.90           N  
ATOM    152  CA  ALA A  33     -13.558  53.413  12.773  1.00 91.90           C  
ATOM    153  C   ALA A  33     -12.043  53.633  12.909  1.00 89.11           C  
ATOM    154  O   ALA A  33     -11.398  54.100  11.970  1.00 83.38           O  
ATOM    155  CB  ALA A  33     -14.366  54.575  13.382  1.00 88.96           C  
ATOM    156  N   ALA A  34     -13.896  49.748  14.975  1.00 87.47           N  
ATOM    157  CA  ALA A  34     -14.087  49.721  13.520  1.00 86.78           C  
ATOM    158  C   ALA A  34     -13.504  50.959  12.879  1.00 89.36           C  
ATOM    159  O   ALA A  34     -12.679  50.860  11.969  1.00 90.48           O  
ATOM    160  CB  ALA A  34     -15.569  49.668  13.136  1.00 86.66           C  
ATOM    161  N   ALA A  35     -13.550  50.695  17.591  1.00 89.35           N  
ATOM    162  CA  ALA A  35     -12.513  49.840  17.010  1.00 90.77           C  
ATOM    163  C   ALA A  35     -12.710  49.495  15.525  1.00 89.57           C  
ATOM    164  O   ALA A  35     -11.787  49.001  14.888  1.00 90.71           O  
ATOM    165  CB  ALA A  35     -12.382  48.551  17.827  1.00 92.76           C  
ATOM    166  N   ALA A  36     -16.091  51.946  17.865  1.00 84.15           N  
ATOM    167  CA  ALA A  36     -14.839  52.725  17.945  1.00 83.59           C  
ATOM    168  C   ALA A  36     -13.671  51.977  17.276  1.00 85.46           C  
ATOM    169  O   ALA A  36     -12.921  52.545  16.480  1.00 82.35           O  
ATOM    170  CB  ALA A  36     -14.481  53.094  19.392  1.00 78.68           C  
TER     171      ALA A  36                                                      
ENDMDL                                                                          
END"""

ENSEMBLE_POLYALA_WITH_FLAGS = """REMARK   2                                                                      
REMARK   2 RESOLUTION.    2.40 ANGSTROMS.                                       
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : REFMAC 5.7.0029                                      
REMARK   3   AUTHORS     : MURSHUDOV,VAGIN,DODSON                               
REMARK   3                                                                      
REMARK   3    REFINEMENT TARGET : MAXIMUM LIKELIHOOD                            
REMARK   3                                                                      
REMARK   3  DATA USED IN REFINEMENT.                                            
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.40                           
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : 44.91                          
REMARK   3   DATA CUTOFF            (SIGMA(F)) : NULL                           
REMARK   3   COMPLETENESS FOR RANGE        (%) : 79.4                           
REMARK   3   NUMBER OF REFLECTIONS             : 9029                           
REMARK   3                                                                      
REMARK   3  FIT TO DATA USED IN REFINEMENT.                                     
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT                      
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM                          
REMARK   3   R VALUE     (WORKING + TEST SET) : 0.198                           
REMARK   3   R VALUE            (WORKING SET) : 0.196                           
REMARK   3   FREE R VALUE                     : 0.248                           
REMARK   3   FREE R VALUE TEST SET SIZE   (%) : 5.000                           
REMARK   3   FREE R VALUE TEST SET COUNT      : 478                             
REMARK   3                                                                      
REMARK   3  FIT IN THE HIGHEST RESOLUTION BIN.                                  
REMARK   3   TOTAL NUMBER OF BINS USED           : 20                           
REMARK   3   BIN RESOLUTION RANGE HIGH       (A) : 2.40                         
REMARK   3   BIN RESOLUTION RANGE LOW        (A) : 2.46                         
REMARK   3   REFLECTION IN BIN     (WORKING SET) : 389                          
REMARK   3   BIN COMPLETENESS (WORKING+TEST) (%) : 47.85                        
REMARK   3   BIN R VALUE           (WORKING SET) : 0.2590                       
REMARK   3   BIN FREE R VALUE SET COUNT          : 23                           
REMARK   3   BIN FREE R VALUE                    : 0.2920                       
REMARK   3                                                                      
REMARK   3  NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.                    
REMARK   3   PROTEIN ATOMS            : 1451                                    
REMARK   3   NUCLEIC ACID ATOMS       : 0                                       
REMARK   3   HETEROGEN ATOMS          : 0                                       
REMARK   3   SOLVENT ATOMS            : 50                                      
REMARK   3                                                                      
REMARK   3  B VALUES.                                                           
REMARK   3   FROM WILSON PLOT           (A**2) : NULL                           
REMARK   3   MEAN B VALUE      (OVERALL, A**2) : 47.53                          
REMARK   3   OVERALL ANISOTROPIC B VALUE.                                       
REMARK   3    B11 (A**2) : -1.71000                                             
REMARK   3    B22 (A**2) : -1.71000                                             
REMARK   3    B33 (A**2) : 5.53000                                              
REMARK   3    B12 (A**2) : -1.71000                                             
REMARK   3    B13 (A**2) : -0.00000                                             
REMARK   3    B23 (A**2) : 0.00000                                              
REMARK   3                                                                      
REMARK   3  ESTIMATED OVERALL COORDINATE ERROR.                                 
REMARK   3   ESU BASED ON R VALUE                            (A): 0.364         
REMARK   3   ESU BASED ON FREE R VALUE                       (A): 0.266         
REMARK   3   ESU BASED ON MAXIMUM LIKELIHOOD                 (A): 0.155         
REMARK   3   ESU FOR B VALUES BASED ON MAXIMUM LIKELIHOOD (A**2): 6.591         
REMARK   3                                                                      
REMARK   3 CORRELATION COEFFICIENTS.                                            
REMARK   3   CORRELATION COEFFICIENT FO-FC      : 0.937                         
REMARK   3   CORRELATION COEFFICIENT FO-FC FREE : 0.908                         
REMARK   3                                                                      
REMARK   3  RMS DEVIATIONS FROM IDEAL VALUES        COUNT    RMS    WEIGHT      
REMARK   3   BOND LENGTHS REFINED ATOMS        (A):  1501 ; 0.018 ; 0.019       
REMARK   3   BOND LENGTHS OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):  2043 ; 1.965 ; 1.924       
REMARK   3   BOND ANGLES OTHERS          (DEGREES):  NULL ;  NULL ;  NULL       
REMARK   3   TORSION ANGLES, PERIOD 1    (DEGREES):   181 ; 7.098 ; 5.000       
REMARK   3   TORSION ANGLES, PERIOD 2    (DEGREES):    59 ;30.230 ;22.203       
REMARK   3   TORSION ANGLES, PERIOD 3    (DEGREES):   233 ;18.954 ;15.000       
REMARK   3   TORSION ANGLES, PERIOD 4    (DEGREES):     6 ;16.340 ;15.000       
REMARK   3   CHIRAL-CENTER RESTRAINTS       (A**3):   219 ; 0.126 ; 0.200       
REMARK   3   GENERAL PLANES REFINED ATOMS      (A):  1127 ; 0.009 ; 0.021       
REMARK   3   GENERAL PLANES OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED CONTACTS OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   NON-BONDED TORSION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) REFINED ATOMS      (A):  NULL ;  NULL ;  NULL       
REMARK   3   H-BOND (X...Y) OTHERS             (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION REFINED ATOMS (A):  NULL ;  NULL ;  NULL       
REMARK   3   POTENTIAL METAL-ION OTHERS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW REFINED ATOMS        (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY VDW OTHERS               (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND REFINED ATOMS     (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY H-BOND OTHERS            (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION REFINED ATOMS  (A):  NULL ;  NULL ;  NULL       
REMARK   3   SYMMETRY METAL-ION OTHERS         (A):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  ISOTROPIC THERMAL FACTOR RESTRAINTS.     COUNT   RMS    WEIGHT      
REMARK   3   MAIN-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   MAIN-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND REFINED ATOMS  (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN BOND OTHER ATOMS    (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE REFINED ATOMS (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SIDE-CHAIN ANGLE OTHER ATOMS   (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B REFINED ATOMS     (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   LONG RANGE B OTHER ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3 ANISOTROPIC THERMAL FACTOR RESTRAINTS.    COUNT   RMS   WEIGHT       
REMARK   3   RIGID-BOND RESTRAINTS          (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; FREE ATOMS         (A**2):  NULL ;  NULL ;  NULL       
REMARK   3   SPHERICITY; BONDED ATOMS       (A**2):  NULL ;  NULL ;  NULL       
REMARK   3                                                                      
REMARK   3  NCS RESTRAINTS STATISTICS                                           
REMARK   3   NUMBER OF DIFFERENT NCS GROUPS : NULL                              
REMARK   3                                                                      
REMARK   3  TLS DETAILS                                                         
REMARK   3   NUMBER OF TLS GROUPS  : NULL                                       
REMARK   3                                                                      
REMARK   3  BULK SOLVENT MODELLING.                                             
REMARK   3   METHOD USED : MASK                                                 
REMARK   3   PARAMETERS FOR MASK CALCULATION                                    
REMARK   3   VDW PROBE RADIUS   : 1.20                                          
REMARK   3   ION PROBE RADIUS   : 0.80                                          
REMARK   3   SHRINKAGE RADIUS   : 0.80                                          
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: NULL                                      
REMARK   4                                                                      
REMARK   4 4NJN COMPLIES WITH FORMAT V. 3.30, 13-JUL-11                             
CRYST1  110.680  110.680  127.710  90.00  90.00 120.00 H 3 2        18          
MODEL        1                                                                  
ATOM      1  N   ALA A   1       6.808  61.263  13.432  1.00 72.61           N  
ANISOU    1  N   ALA A   1     9095   7963  10530    149   1518  -1669       N  
ATOM      2  CA  ALA A   1       7.802  60.901  14.433  1.00 75.34           C  
ANISOU    2  CA  ALA A   1     9317   8124  11185    233   1511  -1496       C  
ATOM      3  C   ALA A   1       8.710  59.751  13.957  1.00 83.54           C  
ANISOU    3  C   ALA A   1    10223   8997  12519    336   1668  -1701       C  
ATOM      4  O   ALA A   1       9.847  59.594  14.417  1.00 89.53           O  
ANISOU    4  O   ALA A   1    10845   9666  13503    427   1706  -1596       O  
ATOM      5  CB  ALA A   1       7.145  60.691  15.789  1.00 75.55           C  
ANISOU    5  CB  ALA A   1     9364   7997  11342    213   1344  -1273       C  
ATOM      6  N   ALA A   2       5.296  61.322  10.847  1.00 71.18           N  
ANISOU    6  N   ALA A   2     9096   8147   9801     20   1606  -2181       N  
ATOM      7  CA  ALA A   2       6.152  62.390  11.390  1.00 71.60           C  
ANISOU    7  CA  ALA A   2     9122   8276   9806     47   1598  -1902       C  
ATOM      8  C   ALA A   2       7.216  61.908  12.350  1.00 76.14           C  
ANISOU    8  C   ALA A   2     9562   8635  10730    135   1627  -1774       C  
ATOM      9  O   ALA A   2       8.425  62.127  12.148  1.00 78.15           O  
ANISOU    9  O   ALA A   2     9721   8935  11036    188   1743  -1750       O  
ATOM     10  CB  ALA A   2       5.322  63.389  12.190  1.00 81.86           C  
ANISOU   10  CB  ALA A   2    10508   9630  10963     -5   1414  -1647       C  
ATOM     11  N   ALA A   3       4.070  58.691  10.381  1.00 78.09           N  
ANISOU   11  N   ALA A   3     9947   8649  11074    -14   1625  -2713       N  
ATOM     12  CA  ALA A   3       4.632  59.476   9.274  1.00 77.89           C  
ANISOU   12  CA  ALA A   3     9952   8925  10717     -7   1737  -2805       C  
ATOM     13  C   ALA A   3       5.589  60.662   9.714  1.00 79.72           C  
ANISOU   13  C   ALA A   3    10158   9278  10852     34   1755  -2495       C  
ATOM     14  O   ALA A   3       6.580  60.956   9.019  1.00 79.47           O  
ANISOU   14  O   ALA A   3    10082   9389  10722     72   1910  -2556       O  
ATOM     15  CB  ALA A   3       3.483  59.989   8.383  1.00 78.46           C  
ANISOU   15  CB  ALA A   3    10144   9263  10404   -112   1654  -2920       C  
ATOM     16  N   ALA A   4       3.150  58.127  13.054  1.00 73.86           N  
ANISOU   16  N   ALA A   4     9370   7685  11008    -44   1372  -2230       N  
ATOM     17  CA  ALA A   4       4.164  57.312  12.379  1.00 79.41           C  
ANISOU   17  CA  ALA A   4     9998   8267  11906     42   1533  -2465       C  
ATOM     18  C   ALA A   4       4.849  58.044  11.235  1.00 79.34           C  
ANISOU   18  C   ALA A   4    10005   8539  11602     70   1658  -2607       C  
ATOM     19  O   ALA A   4       6.064  58.036  11.132  1.00 82.59           O  
ANISOU   19  O   ALA A   4    10329   8935  12116    166   1782  -2605       O  
ATOM     20  CB  ALA A   4       3.529  56.054  11.792  1.00 90.18           C  
ANISOU   20  CB  ALA A   4    11362   9448  13451     -3   1563  -2793       C  
ATOM     21  N   ALA A   5       1.211  60.152  13.195  1.00 76.46           N  
ANISOU   21  N   ALA A   5     9861   8451  10738   -188   1133  -1999       N  
ATOM     22  CA  ALA A   5       2.272  60.073  14.171  1.00 72.86           C  
ANISOU   22  CA  ALA A   5     9339   7838  10506   -107   1163  -1788       C  
ATOM     23  C   ALA A   5       3.430  59.310  13.580  1.00 72.45           C  
ANISOU   23  C   ALA A   5     9214   7668  10645    -28   1316  -1965       C  
ATOM     24  O   ALA A   5       4.556  59.785  13.621  1.00 73.23           O  
ANISOU   24  O   ALA A   5     9266   7815  10740     44   1382  -1871       O  
ATOM     25  CB  ALA A   5       1.830  59.378  15.430  1.00 75.99           C  
ANISOU   25  CB  ALA A   5     9696   8009  11167   -128   1083  -1618       C  
ATOM     26  N   ALA A   6      -0.524  59.587  11.161  1.00 66.38           N  
ANISOU   26  N   ALA A   6     8660   7395   9166   -337   1099  -2523       N  
ATOM     27  CA  ALA A   6       0.233  60.855  11.098  1.00 63.60           C  
ANISOU   27  CA  ALA A   6     8346   7236   8581   -274   1122  -2331       C  
ATOM     28  C   ALA A   6       1.411  60.765  12.032  1.00 67.49           C  
ANISOU   28  C   ALA A   6     8775   7548   9320   -188   1188  -2138       C  
ATOM     29  O   ALA A   6       2.505  61.188  11.690  1.00 69.06           O  
ANISOU   29  O   ALA A   6     8959   7820   9461   -127   1289  -2114       O  
ATOM     30  CB  ALA A   6      -0.614  62.070  11.372  1.00 56.93           C  
ANISOU   30  CB  ALA A   6     7562   6584   7483   -304    980  -2128       C  
ATOM     31  N   ALA A   7      -1.601  57.282  11.899  1.00 61.88           N  
ANISOU   31  N   ALA A   7     7992   6336   9182   -453   1058  -2717       N  
ATOM     32  CA  ALA A   7      -0.766  57.188  10.752  1.00 67.84           C  
ANISOU   32  CA  ALA A   7     8764   7177   9832   -397   1183  -2972       C  
ATOM     33  C   ALA A   7       0.013  58.494  10.624  1.00 70.11           C  
ANISOU   33  C   ALA A   7     9092   7710   9835   -319   1216  -2784       C  
ATOM     34  O   ALA A   7       1.136  58.503  10.080  1.00 69.35           O  
ANISOU   34  O   ALA A   7     8979   7645   9722   -241   1355  -2874       O  
ATOM     35  CB  ALA A   7      -1.441  56.665   9.450  1.00 73.01           C  
ANISOU   35  CB  ALA A   7     9448   7939  10350   -483   1185  -3372       C  
ATOM     36  N   ALA A   8      -2.856  58.098  14.278  1.00 60.88           N  
ANISOU   36  N   ALA A   8     7836   6188   9105   -513    856  -2123       N  
ATOM     37  CA  ALA A   8      -1.791  57.132  14.351  1.00 59.72           C  
ANISOU   37  CA  ALA A   8     7651   5777   9263   -452    963  -2180       C  
ATOM     38  C   ALA A   8      -0.982  57.170  13.066  1.00 60.98           C  
ANISOU   38  C   ALA A   8     7835   6031   9301   -391   1071  -2449       C  
ATOM     39  O   ALA A   8       0.223  57.129  13.131  1.00 62.29           O  
ANISOU   39  O   ALA A   8     7974   6123   9569   -287   1166  -2409       O  
ATOM     40  CB  ALA A   8      -2.352  55.726  14.631  1.00 62.50           C  
ANISOU   40  CB  ALA A   8     7947   5828   9972   -539    960  -2276       C  
ATOM     41  N   ALA A   9      -4.846  59.990  13.395  1.00 58.40           N  
ANISOU   41  N   ALA A   9     7589   6436   8161   -595    652  -2149       N  
ATOM     42  CA  ALA A   9      -3.712  60.360  14.217  1.00 60.43           C  
ANISOU   42  CA  ALA A   9     7859   6604   8495   -498    713  -1918       C  
ATOM     43  C   ALA A   9      -2.570  59.381  14.182  1.00 58.45           C  
ANISOU   43  C   ALA A   9     7580   6117   8509   -453    834  -1997       C  
ATOM     44  O   ALA A   9      -1.430  59.784  14.089  1.00 59.36           O  
ANISOU   44  O   ALA A   9     7710   6250   8590   -363    905  -1939       O  
ATOM     45  CB  ALA A   9      -4.092  60.532  15.675  1.00 58.59           C  
ANISOU   45  CB  ALA A   9     7590   6300   8370   -505    657  -1640       C  
ATOM     46  N   ALA A  10      -6.665  58.234  12.008  1.00 64.10           N  
ANISOU   46  N   ALA A  10     8232   7113   9010   -830    573  -2697       N  
ATOM     47  CA  ALA A  10      -5.894  59.264  11.295  1.00 66.28           C  
ANISOU   47  CA  ALA A  10     8594   7613   8975   -726    604  -2683       C  
ATOM     48  C   ALA A  10      -4.683  59.763  12.098  1.00 63.40           C  
ANISOU   48  C   ALA A  10     8255   7161   8671   -608    688  -2424       C  
ATOM     49  O   ALA A  10      -3.609  59.922  11.539  1.00 68.98           O  
ANISOU   49  O   ALA A  10     9001   7896   9309   -536    786  -2482       O  
ATOM     50  CB  ALA A  10      -6.774  60.431  10.945  1.00 62.41           C  
ANISOU   50  CB  ALA A  10     8127   7431   8153   -729    472  -2610       C  
ATOM     51  N   ALA A  11      -7.523  56.981  14.399  1.00 63.45           N  
ANISOU   51  N   ALA A  11     8001   6600   9506   -957    560  -2345       N  
ATOM     52  CA  ALA A  11      -6.926  56.234  13.306  1.00 64.48           C  
ANISOU   52  CA  ALA A  11     8161   6628   9707   -955    625  -2661       C  
ATOM     53  C   ALA A  11      -6.072  57.139  12.464  1.00 63.66           C  
ANISOU   53  C   ALA A  11     8144   6743   9299   -835    664  -2722       C  
ATOM     54  O   ALA A  11      -4.886  56.872  12.245  1.00 60.49           O  
ANISOU   54  O   ALA A  11     7771   6229   8982   -747    777  -2775       O  
ATOM     55  N   ALA A  12      -8.666  58.912  15.975  1.00 61.98           N  
ANISOU   55  N   ALA A  12     7779   6763   9007   -914    436  -1866       N  
ATOM     56  CA  ALA A  12      -7.453  58.270  16.446  1.00 60.60           C  
ANISOU   56  CA  ALA A  12     7634   6335   9053   -862    535  -1781       C  
ATOM     57  C   ALA A  12      -6.762  57.549  15.313  1.00 64.32           C  
ANISOU   57  C   ALA A  12     8139   6692   9606   -853    596  -2064       C  
ATOM     58  O   ALA A  12      -5.548  57.496  15.275  1.00 65.65           O  
ANISOU   58  O   ALA A  12     8348   6762   9831   -754    676  -2042       O  
ATOM     59  CB  ALA A  12      -7.713  57.288  17.567  1.00 60.13           C  
ANISOU   59  CB  ALA A  12     7504   6036   9306   -947    562  -1617       C  
ATOM     60  N   ALA A  13     -10.929  59.546  14.435  1.00 68.55           N  
ANISOU   60  N   ALA A  13     8519   8026   9501  -1036    226  -2184       N  
ATOM     61  CA  ALA A  13      -9.869  60.510  14.575  1.00 63.92           C  
ANISOU   61  CA  ALA A  13     8040   7496   8748   -885    267  -2029       C  
ATOM     62  C   ALA A  13      -8.582  59.839  15.040  1.00 63.57           C  
ANISOU   62  C   ALA A  13     8042   7192   8919   -846    392  -1965       C  
ATOM     63  O   ALA A  13      -7.516  60.131  14.527  1.00 72.56           O  
ANISOU   63  O   ALA A  13     9262   8339   9969   -758    447  -2000       O  
ATOM     64  CB  ALA A  13     -10.261  61.598  15.526  1.00 62.67           C  
ANISOU   64  CB  ALA A  13     7865   7469   8476   -815    225  -1776       C  
ATOM     65  N   ALA A  14     -12.378  57.236  14.902  1.00 66.70           N  
ANISOU   65  N   ALA A  14     8073   7456   9813  -1356    227  -2355       N  
ATOM     66  CA  ALA A  14     -11.949  57.558  13.541  1.00 67.37           C  
ANISOU   66  CA  ALA A  14     8246   7680   9670  -1296    192  -2604       C  
ATOM     67  C   ALA A  14     -10.808  58.532  13.580  1.00 65.65           C  
ANISOU   67  C   ALA A  14     8152   7540   9251  -1118    245  -2455       C  
ATOM     68  O   ALA A  14      -9.829  58.360  12.875  1.00 64.93           O  
ANISOU   68  O   ALA A  14     8148   7399   9123  -1062    310  -2590       O  
ATOM     69  CB  ALA A  14     -13.052  58.190  12.723  1.00 66.46           C  
ANISOU   69  CB  ALA A  14     8075   7875   9300  -1332     45  -2738       C  
ATOM     70  N   ALA A  15     -12.829  57.552  17.677  1.00 70.33           N  
ANISOU   70  N   ALA A  15     8426   7852  10442  -1361    311  -1747       N  
ATOM     71  CA  ALA A  15     -12.001  56.473  17.189  1.00 69.48           C  
ANISOU   71  CA  ALA A  15     8374   7462  10560  -1392    366  -1894       C  
ATOM     72  C   ALA A  15     -11.494  56.785  15.792  1.00 67.75           C  
ANISOU   72  C   ALA A  15     8244   7344  10151  -1318    333  -2170       C  
ATOM     73  O   ALA A  15     -10.324  56.588  15.510  1.00 67.50           O  
ANISOU   73  O   ALA A  15     8306   7181  10160  -1230    403  -2214       O  
ATOM     74  N   ALA A  16     -14.520  59.746  17.510  1.00 64.73           N  
ANISOU   74  N   ALA A  16     7610   7736   9247  -1259    144  -1711       N  
ATOM     75  CA  ALA A  16     -13.268  59.872  18.256  1.00 60.88           C  
ANISOU   75  CA  ALA A  16     7243   7104   8784  -1159    239  -1519       C  
ATOM     76  C   ALA A  16     -12.305  58.751  17.863  1.00 66.60           C  
ANISOU   76  C   ALA A  16     8035   7547   9720  -1201    303  -1621       C  
ATOM     77  O   ALA A  16     -11.112  58.968  17.689  1.00 64.53           O  
ANISOU   77  O   ALA A  16     7889   7216   9412  -1090    345  -1601       O  
ATOM     78  CB  ALA A  16     -13.497  59.841  19.762  1.00 58.66           C  
ANISOU   78  CB  ALA A  16     6902   6796   8590  -1182    308  -1249       C  
ATOM     79  N   ALA A  17      -0.741  56.623  -1.910  1.00 92.22           N  
ANISOU   79  N   ALA A  17    12323  12754   9960   -712   1738  -6044       N  
ATOM     80  CA  ALA A  17      -1.060  57.039  -3.284  1.00 97.80           C  
ANISOU   80  CA  ALA A  17    13124  13951  10083   -788   1712  -6239       C  
ATOM     81  C   ALA A  17      -0.238  58.256  -3.688  1.00 96.58           C  
ANISOU   81  C   ALA A  17    13022  14135   9538   -734   1831  -5932       C  
ATOM     82  O   ALA A  17       0.400  58.251  -4.743  1.00100.43           O  
ANISOU   82  O   ALA A  17    13544  14922   9690   -729   2010  -6157       O  
ATOM     83  CB  ALA A  17      -2.585  57.224  -3.546  1.00 97.39           C  
ANISOU   83  CB  ALA A  17    13115  14073   9813   -923   1410  -6256       C  
ATOM     84  N   ALA A  18       0.124  56.418   0.683  1.00 90.78           N  
ANISOU   84  N   ALA A  18    11988  11815  10687   -550   1784  -5476       N  
ATOM     85  CA  ALA A  18       0.636  55.459  -0.268  1.00 89.62           C  
ANISOU   85  CA  ALA A  18    11827  11678  10545   -533   1962  -5969       C  
ATOM     86  C   ALA A  18       0.389  55.958  -1.674  1.00 92.73           C  
ANISOU   86  C   ALA A  18    12316  12569  10345   -605   1968  -6175       C  
ATOM     87  O   ALA A  18       1.219  55.731  -2.548  1.00100.14           O  
ANISOU   87  O   ALA A  18    13256  13678  11116   -559   2182  -6451       O  
ATOM     88  CB  ALA A  18       0.020  54.059  -0.079  1.00 91.24           C  
ANISOU   88  CB  ALA A  18    11989  11521  11157   -589   1904  -6320       C  
ATOM     89  N   ALA A  19      -1.443  58.421   1.789  1.00 83.59           N  
ANISOU   89  N   ALA A  19    11148  11110   9501   -629   1388  -4720       N  
ATOM     90  CA  ALA A  19      -0.010  58.556   1.757  1.00 83.29           C  
ANISOU   90  CA  ALA A  19    11087  11051   9509   -521   1614  -4669       C  
ATOM     91  C   ALA A  19       0.613  57.651   0.742  1.00 87.40           C  
ANISOU   91  C   ALA A  19    11593  11616   9996   -506   1808  -5122       C  
ATOM     92  O   ALA A  19       1.522  58.047   0.027  1.00 91.61           O  
ANISOU   92  O   ALA A  19    12141  12383  10283   -460   1983  -5160       O  
ATOM     93  N   ALA A  20      -4.028  57.264   1.493  1.00 82.93           N  
ANISOU   93  N   ALA A  20    11039  10986   9483   -869   1030  -5124       N  
ATOM     94  CA  ALA A  20      -3.666  58.563   0.892  1.00 82.67           C  
ANISOU   94  CA  ALA A  20    11088  11335   8987   -821   1041  -4909       C  
ATOM     95  C   ALA A  20      -2.168  58.728   0.737  1.00 84.53           C  
ANISOU   95  C   ALA A  20    11336  11571   9211   -712   1288  -4874       C  
ATOM     96  O   ALA A  20      -1.677  59.143  -0.306  1.00 93.03           O  
ANISOU   96  O   ALA A  20    12472  12966   9906   -705   1385  -4969       O  
ATOM     97  CB  ALA A  20      -4.188  59.702   1.728  1.00 79.90           C  
ANISOU   97  CB  ALA A  20    10741  11020   8596   -799    885  -4451       C  
ATOM     98  N   ALA A  21      -3.959  54.929   3.149  1.00 84.31           N  
ANISOU   98  N   ALA A  21    11053  10268  10712   -888   1110  -5307       N  
ATOM     99  CA  ALA A  21      -4.144  54.839   1.692  1.00 87.02           C  
ANISOU   99  CA  ALA A  21    11454  10928  10679   -948   1113  -5701       C  
ATOM    100  C   ALA A  21      -3.663  56.089   0.966  1.00 86.80           C  
ANISOU  100  C   ALA A  21    11509  11332  10136   -886   1153  -5548       C  
ATOM    101  O   ALA A  21      -2.986  55.978  -0.066  1.00 89.18           O  
ANISOU  101  O   ALA A  21    11855  11833  10196   -859   1299  -5818       O  
ATOM    102  CB  ALA A  21      -5.602  54.562   1.317  1.00 87.25           C  
ANISOU  102  CB  ALA A  21    11469  11063  10616  -1114    893  -5889       C  
ATOM    103  N   ALA A  22      -3.535  56.148   5.801  1.00 84.81           N  
ANISOU  103  N   ALA A  22    11058  10047  11117   -750   1063  -4375       N  
ATOM    104  CA  ALA A  22      -2.598  55.218   5.170  1.00 85.75           C  
ANISOU  104  CA  ALA A  22    11166  10027  11385   -700   1242  -4709       C  
ATOM    105  C   ALA A  22      -2.760  55.231   3.639  1.00 88.21           C  
ANISOU  105  C   ALA A  22    11540  10663  11311   -751   1264  -5093       C  
ATOM    106  O   ALA A  22      -1.810  55.551   2.908  1.00 88.95           O  
ANISOU  106  O   ALA A  22    11668  10944  11183   -672   1417  -5184       O  
ATOM    107  CB  ALA A  22      -2.820  53.811   5.697  1.00 85.65           C  
ANISOU  107  CB  ALA A  22    11084   9586  11874   -746   1250  -4875       C  
ATOM    108  N   ALA A  23      -6.157  57.136   5.909  1.00 79.85           N  
ANISOU  108  N   ALA A  23    10417   9768  10154   -941    676  -4199       N  
ATOM    109  CA  ALA A  23      -5.050  58.046   5.818  1.00 83.06           C  
ANISOU  109  CA  ALA A  23    10889  10305  10365   -812    784  -4009       C  
ATOM    110  C   ALA A  23      -3.867  57.306   5.229  1.00 85.48           C  
ANISOU  110  C   ALA A  23    11206  10506  10764   -756    983  -4270       C  
ATOM    111  O   ALA A  23      -3.281  57.777   4.258  1.00 83.50           O  
ANISOU  111  O   ALA A  23    11016  10514  10194   -717   1065  -4372       O  
ATOM    112  CB  ALA A  23      -4.772  58.617   7.180  1.00 81.77           C  
ANISOU  112  CB  ALA A  23    10700   9991  10376   -740    781  -3602       C  
ATOM    113  N   ALA A  24      -7.769  54.771   5.827  1.00 78.19           N  
ANISOU  113  N   ALA A  24    10068   9163  10478  -1217    559  -4733       N  
ATOM    114  CA  ALA A  24      -8.046  55.947   4.990  1.00 81.68           C  
ANISOU  114  CA  ALA A  24    10571  10049  10412  -1191    468  -4698       C  
ATOM    115  C   ALA A  24      -6.839  56.880   4.827  1.00 81.93           C  
ANISOU  115  C   ALA A  24    10692  10228  10207  -1035    595  -4511       C  
ATOM    116  O   ALA A  24      -6.513  57.366   3.733  1.00 89.78           O  
ANISOU  116  O   ALA A  24    11761  11522  10827  -1010    620  -4642       O  
ATOM    117  CB  ALA A  24      -9.211  56.733   5.574  1.00 81.22           C  
ANISOU  117  CB  ALA A  24    10456  10136  10264  -1234    285  -4421       C  
ATOM    118  N   ALA A  25      -6.824  53.011   7.753  1.00 75.58           N  
ANISOU  118  N   ALA A  25     9643   7996  11077  -1182    745  -4548       N  
ATOM    119  CA  ALA A  25      -6.776  52.613   6.349  1.00 78.65           C  
ANISOU  119  CA  ALA A  25    10072   8518  11291  -1219    767  -5010       C  
ATOM    120  C   ALA A  25      -6.915  53.834   5.435  1.00 79.93           C  
ANISOU  120  C   ALA A  25    10305   9169  10894  -1187    707  -5011       C  
ATOM    121  O   ALA A  25      -6.198  53.953   4.425  1.00 80.99           O  
ANISOU  121  O   ALA A  25    10506   9481  10783  -1126    804  -5239       O  
ATOM    122  CB  ALA A  25      -7.833  51.534   6.048  1.00 79.68           C  
ANISOU  122  CB  ALA A  25    10141   8491  11642  -1409    668  -5329       C  
ATOM    123  N   ALA A  26      -7.484  54.579   9.963  1.00 74.92           N  
ANISOU  123  N   ALA A  26     9509   8016  10940  -1136    621  -3715       N  
ATOM    124  CA  ALA A  26      -6.092  54.199   9.754  1.00 71.92           C  
ANISOU  124  CA  ALA A  26     9173   7471  10680  -1015    777  -3797       C  
ATOM    125  C   ALA A  26      -5.907  53.803   8.270  1.00 70.77           C  
ANISOU  125  C   ALA A  26     9070   7439  10380  -1034    821  -4249       C  
ATOM    126  O   ALA A  26      -4.995  54.258   7.598  1.00 65.93           O  
ANISOU  126  O   ALA A  26     8515   6986   9547   -927    919  -4325       O  
ATOM    127  CB  ALA A  26      -5.649  53.118  10.737  1.00 68.53           C  
ANISOU  127  CB  ALA A  26     8689   6606  10742  -1009    849  -3690       C  
ATOM    128  N   ALA A  27     -10.265  54.880   9.568  1.00 76.52           N  
ANISOU  128  N   ALA A  27     9582   8571  10921  -1403    317  -3833       N  
ATOM    129  CA  ALA A  27      -9.411  56.032   9.567  1.00 73.50           C  
ANISOU  129  CA  ALA A  27     9291   8376  10259  -1234    364  -3631       C  
ATOM    130  C   ALA A  27      -7.957  55.661   9.340  1.00 75.75           C  
ANISOU  130  C   ALA A  27     9650   8495  10636  -1125    529  -3709       C  
ATOM    131  O   ALA A  27      -7.258  56.344   8.580  1.00 74.64           O  
ANISOU  131  O   ALA A  27     9590   8566  10203  -1030    575  -3760       O  
ATOM    132  N   ALA A  28     -10.888  52.210   9.995  1.00 79.48           N  
ANISOU  132  N   ALA A  28     9818   8278  12100  -1678    353  -4147       N  
ATOM    133  CA  ALA A  28     -11.295  52.880   8.749  1.00 81.36           C  
ANISOU  133  CA  ALA A  28    10090   8905  11917  -1683    253  -4393       C  
ATOM    134  C   ALA A  28     -10.430  54.097   8.506  1.00 78.92           C  
ANISOU  134  C   ALA A  28     9887   8857  11238  -1493    301  -4233       C  
ATOM    135  O   ALA A  28      -9.986  54.334   7.388  1.00 80.63           O  
ANISOU  135  O   ALA A  28    10183   9270  11180  -1445    318  -4468       O  
ATOM    136  CB  ALA A  28     -12.753  53.363   8.798  1.00 81.84           C  
ANISOU  136  CB  ALA A  28    10043   9218  11834  -1803     77  -4332       C  
ATOM    137  N   ALA A  29      -9.921  51.628  12.646  1.00 71.50           N  
ANISOU  137  N   ALA A  29     8772   6704  11688  -1591    519  -3453       N  
ATOM    138  CA  ALA A  29      -9.333  50.983  11.475  1.00 76.71           C  
ANISOU  138  CA  ALA A  29     9485   7265  12394  -1572    569  -3863       C  
ATOM    139  C   ALA A  29      -9.682  51.674  10.130  1.00 78.29           C  
ANISOU  139  C   ALA A  29     9727   7858  12159  -1574    494  -4163       C  
ATOM    140  O   ALA A  29      -8.854  51.733   9.230  1.00 84.74           O  
ANISOU  140  O   ALA A  29    10623   8746  12825  -1478    565  -4395       O  
ATOM    141  CB  ALA A  29      -9.735  49.504  11.425  1.00 80.04           C  
ANISOU  141  CB  ALA A  29     9844   7299  13266  -1734    571  -4087       C  
ATOM    142  N   ALA A  30     -11.893  53.345  13.999  1.00 74.27           N  
ANISOU  142  N   ALA A  30     8986   7559  11672  -1673    352  -2939       N  
ATOM    143  CA  ALA A  30     -10.435  53.512  14.126  1.00 71.09           C  
ANISOU  143  CA  ALA A  30     8697   7047  11264  -1497    457  -2848       C  
ATOM    144  C   ALA A  30      -9.670  52.899  12.933  1.00 71.27           C  
ANISOU  144  C   ALA A  30     8787   6960  11331  -1463    508  -3215       C  
ATOM    145  O   ALA A  30      -8.873  53.580  12.288  1.00 68.66           O  
ANISOU  145  O   ALA A  30     8543   6802  10742  -1328    544  -3285       O  
ATOM    146  CB  ALA A  30      -9.948  52.873  15.427  1.00 69.28           C  
ANISOU  146  CB  ALA A  30     8448   6504  11369  -1492    537  -2547       C  
ATOM    147  N   ALA A  31     -14.263  51.971  13.343  1.00 79.28           N  
ANISOU  147  N   ALA A  31     9380   8145  12597  -2085    186  -3327       N  
ATOM    148  CA  ALA A  31     -14.060  53.330  12.895  1.00 75.44           C  
ANISOU  148  CA  ALA A  31     8965   8031  11668  -1920    139  -3291       C  
ATOM    149  C   ALA A  31     -12.565  53.670  12.892  1.00 75.27           C  
ANISOU  149  C   ALA A  31     9097   7939  11563  -1721    249  -3212       C  
ATOM    150  O   ALA A  31     -12.038  54.210  11.900  1.00 76.31           O  
ANISOU  150  O   ALA A  31     9321   8256  11417  -1620    239  -3387       O  
ATOM    151  CB  ALA A  31     -14.840  54.279  13.775  1.00 70.59           C  
ANISOU  151  CB  ALA A  31     8265   7650  10906  -1899     95  -2983       C  
ATOM    152  N   ALA A  32     -13.771  49.624  14.718  1.00 81.52           N  
ANISOU  152  N   ALA A  32     9616   7596  13759  -2270    348  -3149       N  
ATOM    153  CA  ALA A  32     -13.901  49.565  13.278  1.00 83.38           C  
ANISOU  153  CA  ALA A  32     9879   7954  13845  -2292    276  -3601       C  
ATOM    154  C   ALA A  32     -13.629  50.972  12.753  1.00 81.26           C  
ANISOU  154  C   ALA A  32     9688   8113  13073  -2112    236  -3592       C  
ATOM    155  O   ALA A  32     -12.785  51.164  11.877  1.00 87.42           O  
ANISOU  155  O   ALA A  32    10580   8952  13683  -1990    266  -3802       O  
ATOM    156  CB  ALA A  32     -15.266  48.953  12.847  1.00 83.45           C  
ANISOU  156  CB  ALA A  32     9743   7974  13989  -2548    161  -3838       C  
TER     157      ALA A  32                                                      
ENDMDL                                                                          
MODEL        2                                                                  
ATOM      1  N   ALA A   1       8.058  59.006  12.665  1.00 80.38           N  
ATOM      2  CA  ALA A   1       8.986  57.987  12.141  1.00 83.37           C  
ATOM      3  C   ALA A   1       9.621  58.375  10.813  1.00 84.33           C  
ATOM      4  O   ALA A   1      10.857  58.412  10.667  1.00 72.14           O  
ATOM      5  N   ALA A   2       6.399  61.381  12.666  1.00 84.79           N  
ATOM      6  CA  ALA A   2       7.486  61.188  13.644  1.00 83.84           C  
ATOM      7  C   ALA A   2       8.506  60.185  13.086  1.00 82.59           C  
ATOM      8  O   ALA A   2       9.690  60.499  13.007  1.00 80.71           O  
ATOM      9  CB  ALA A   2       7.005  60.800  15.065  1.00 81.10           C  
ATOM     10  N   ALA A   3       4.787  61.058  10.259  1.00 98.19           N  
ATOM     11  CA  ALA A   3       5.462  62.335  10.588  1.00 94.68           C  
ATOM     12  C   ALA A   3       6.542  62.260  11.678  1.00 90.57           C  
ATOM     13  O   ALA A   3       7.503  63.029  11.617  1.00 95.24           O  
ATOM     14  CB  ALA A   3       4.433  63.416  10.920  1.00 97.47           C  
ATOM     15  N   ALA A   4       3.940  58.236  10.061  1.00 98.28           N  
ATOM     16  CA  ALA A   4       4.240  59.092   8.895  1.00100.67           C  
ATOM     17  C   ALA A   4       5.056  60.365   9.148  1.00103.54           C  
ATOM     18  O   ALA A   4       5.935  60.702   8.349  1.00112.56           O  
ATOM     19  N   ALA A   5       3.404  57.550  12.915  1.00 85.39           N  
ATOM     20  CA  ALA A   5       4.372  56.881  12.035  1.00 89.96           C  
ATOM     21  C   ALA A   5       4.868  57.794  10.915  1.00 93.64           C  
ATOM     22  O   ALA A   5       6.066  58.073  10.832  1.00 94.85           O  
ATOM     23  CB  ALA A   5       3.792  55.605  11.413  1.00 95.55           C  
ATOM     24  N   ALA A   6       1.500  59.606  13.693  1.00 86.86           N  
ATOM     25  CA  ALA A   6       2.680  59.273  14.521  1.00 87.89           C  
ATOM     26  C   ALA A   6       3.783  58.508  13.760  1.00 87.20           C  
ATOM     27  O   ALA A   6       4.965  58.775  13.962  1.00 84.49           O  
ATOM     28  CB  ALA A   6       2.259  58.462  15.751  1.00 90.68           C  
ATOM     29  N   ALA A   7      -0.503  59.130  11.667  1.00 88.27           N  
ATOM     30  CA  ALA A   7       0.320  60.338  11.651  1.00 86.96           C  
ATOM     31  C   ALA A   7       1.589  60.245  12.522  1.00 86.52           C  
ATOM     32  O   ALA A   7       2.635  60.742  12.115  1.00 86.56           O  
ATOM     33  CB  ALA A   7      -0.513  61.547  12.041  1.00 90.39           C  
ATOM     34  N   ALA A   8      -1.759  56.997  13.061  1.00 88.85           N  
ATOM     35  CA  ALA A   8      -0.936  56.773  11.879  1.00 91.83           C  
ATOM     36  C   ALA A   8       0.028  57.910  11.669  1.00 91.48           C  
ATOM     37  O   ALA A   8       1.228  57.681  11.516  1.00 89.90           O  
ATOM     38  CB  ALA A   8      -1.800  56.611  10.635  1.00 97.68           C  
ATOM     39  N   ALA A   9      -3.104  58.390  15.089  1.00 76.97           N  
ATOM     40  CA  ALA A   9      -2.224  57.257  15.420  1.00 80.80           C  
ATOM     41  C   ALA A   9      -1.249  56.959  14.290  1.00 84.65           C  
ATOM     42  O   ALA A   9      -0.069  56.694  14.533  1.00 80.30           O  
ATOM     43  CB  ALA A   9      -3.021  55.986  15.739  1.00 83.11           C  
ATOM     44  N   ALA A  10      -4.659  60.035  13.347  1.00 78.44           N  
ATOM     45  CA  ALA A  10      -3.647  60.607  14.243  1.00 78.36           C  
ATOM     46  C   ALA A  10      -2.624  59.560  14.670  1.00 76.40           C  
ATOM     47  O   ALA A  10      -1.428  59.821  14.646  1.00 77.82           O  
ATOM     48  CB  ALA A  10      -4.310  61.193  15.496  1.00 80.79           C  
ATOM     49  N   ALA A  11      -6.392  58.056  12.300  1.00 80.31           N  
ATOM     50  CA  ALA A  11      -5.555  58.884  11.424  1.00 78.65           C  
ATOM     51  C   ALA A  11      -4.385  59.539  12.140  1.00 80.24           C  
ATOM     52  O   ALA A  11      -3.268  59.607  11.599  1.00 81.32           O  
ATOM     53  N   ALA A  12      -7.617  57.158  14.794  1.00 79.52           N  
ATOM     54  CA  ALA A  12      -6.969  56.260  13.822  1.00 83.16           C  
ATOM     55  C   ALA A  12      -5.967  56.929  12.882  1.00 83.52           C  
ATOM     56  O   ALA A  12      -4.845  56.424  12.690  1.00 80.28           O  
ATOM     57  N   ALA A  13      -8.829  59.390  16.057  1.00 80.77           N  
ATOM     58  CA  ALA A  13      -7.674  58.730  16.671  1.00 78.56           C  
ATOM     59  C   ALA A  13      -6.912  57.819  15.707  1.00 81.26           C  
ATOM     60  O   ALA A  13      -5.691  57.712  15.806  1.00 87.75           O  
ATOM     61  CB  ALA A  13      -8.096  57.930  17.909  1.00 77.34           C  
ATOM     62  N   ALA A  14     -10.938  59.775  14.113  1.00 85.16           N  
ATOM     63  CA  ALA A  14      -9.992  60.824  14.486  1.00 85.86           C  
ATOM     64  C   ALA A  14      -8.704  60.226  15.027  1.00 87.77           C  
ATOM     65  O   ALA A  14      -7.622  60.494  14.496  1.00 94.35           O  
ATOM     66  N   ALA A  15     -12.188  57.250  14.520  1.00 74.32           N  
ATOM     67  CA  ALA A  15     -11.660  57.640  13.196  1.00 78.10           C  
ATOM     68  C   ALA A  15     -10.614  58.752  13.320  1.00 80.86           C  
ATOM     69  O   ALA A  15      -9.544  58.673  12.717  1.00 80.60           O  
ATOM     70  CB  ALA A  15     -12.823  58.101  12.304  1.00 78.18           C  
ATOM     71  N   ALA A  16     -13.313  57.509  17.131  1.00 67.46           N  
ATOM     72  CA  ALA A  16     -12.142  56.652  16.874  1.00 67.58           C  
ATOM     73  C   ALA A  16     -11.418  56.905  15.551  1.00 69.32           C  
ATOM     74  O   ALA A  16     -10.187  56.777  15.461  1.00 63.60           O  
ATOM     75  N   ALA A  17     -15.707  59.079  16.716  1.00 78.72           N  
ATOM     76  CA  ALA A  17     -14.532  59.651  17.408  1.00 76.19           C  
ATOM     77  C   ALA A  17     -13.211  58.825  17.319  1.00 71.43           C  
ATOM     78  O   ALA A  17     -12.123  59.397  17.439  1.00 70.66           O  
ATOM     79  CB  ALA A  17     -14.877  60.042  18.883  1.00 74.40           C  
ATOM     80  N   ALA A  18      -0.307  59.327  -2.737  1.00 67.33           N  
ATOM     81  CA  ALA A  18       0.496  60.540  -2.776  1.00 67.32           C  
ATOM     82  C   ALA A  18       1.975  60.281  -2.974  1.00 68.67           C  
ATOM     83  O   ALA A  18       2.629  60.972  -3.754  1.00 70.56           O  
ATOM     84  N   ALA A  19      -0.992  56.651  -2.261  1.00 71.55           N  
ATOM     85  CA  ALA A  19      -1.326  57.287  -3.542  1.00 71.39           C  
ATOM     86  C   ALA A  19      -0.524  58.558  -3.803  1.00 70.39           C  
ATOM     87  O   ALA A  19      -0.119  58.819  -4.940  1.00 71.13           O  
ATOM     88  CB  ALA A  19      -2.812  57.621  -3.590  1.00 72.00           C  
ATOM     89  N   ALA A  20      -0.104  56.693   0.468  1.00 86.19           N  
ATOM     90  CA  ALA A  20       0.440  55.756  -0.534  1.00 78.57           C  
ATOM     91  C   ALA A  20       0.252  56.304  -1.934  1.00 76.29           C  
ATOM     92  O   ALA A  20       1.210  56.384  -2.700  1.00 80.88           O  
ATOM     93  CB  ALA A  20      -0.235  54.383  -0.447  1.00 78.63           C  
ATOM     94  N   ALA A  21      -1.842  58.649   1.587  1.00 84.05           N  
ATOM     95  CA  ALA A  21      -0.383  58.829   1.618  1.00 87.23           C  
ATOM     96  C   ALA A  21       0.391  57.921   0.664  1.00 90.58           C  
ATOM     97  O   ALA A  21       1.438  58.325   0.134  1.00 94.29           O  
ATOM     98  N   ALA A  22      -4.232  57.284   1.260  1.00 86.55           N  
ATOM     99  CA  ALA A  22      -4.063  58.580   0.586  1.00 85.23           C  
ATOM    100  C   ALA A  22      -2.606  59.036   0.562  1.00 84.66           C  
ATOM    101  O   ALA A  22      -2.197  59.718  -0.373  1.00 89.02           O  
ATOM    102  CB  ALA A  22      -4.941  59.652   1.219  1.00 85.10           C  
ATOM    103  N   ALA A  23      -3.693  54.944   2.969  1.00 75.30           N  
ATOM    104  CA  ALA A  23      -3.905  54.875   1.529  1.00 79.24           C  
ATOM    105  C   ALA A  23      -3.611  56.194   0.816  1.00 84.91           C  
ATOM    106  O   ALA A  23      -2.818  56.220  -0.125  1.00 88.99           O  
ATOM    107  CB  ALA A  23      -5.328  54.375   1.208  1.00 81.24           C  
ATOM    108  N   ALA A  24      -3.241  55.970   5.751  1.00 78.32           N  
ATOM    109  CA  ALA A  24      -2.388  55.028   4.998  1.00 80.32           C  
ATOM    110  C   ALA A  24      -2.499  55.183   3.488  1.00 78.29           C  
ATOM    111  O   ALA A  24      -1.523  55.538   2.825  1.00 81.62           O  
ATOM    112  CB  ALA A  24      -2.696  53.566   5.386  1.00 80.86           C  
ATOM    113  N   ALA A  25      -5.653  57.311   6.146  1.00 75.64           N  
ATOM    114  CA  ALA A  25      -4.367  58.045   6.269  1.00 75.50           C  
ATOM    115  C   ALA A  25      -3.355  57.264   5.442  1.00 77.16           C  
ATOM    116  O   ALA A  25      -2.743  57.812   4.521  1.00 79.43           O  
ATOM    117  CB  ALA A  25      -3.895  58.185   7.747  1.00 74.94           C  
ATOM    118  N   ALA A  26      -7.439  55.056   5.824  1.00 65.86           N  
ATOM    119  CA  ALA A  26      -7.523  56.225   4.945  1.00 73.62           C  
ATOM    120  C   ALA A  26      -6.256  57.116   4.967  1.00 74.88           C  
ATOM    121  O   ALA A  26      -5.843  57.617   3.914  1.00 75.74           O  
ATOM    122  CB  ALA A  26      -8.766  57.065   5.264  1.00 77.48           C  
ATOM    123  N   ALA A  27      -6.580  53.403   7.935  1.00 65.93           N  
ATOM    124  CA  ALA A  27      -6.464  52.939   6.562  1.00 64.22           C  
ATOM    125  C   ALA A  27      -6.558  54.105   5.577  1.00 65.79           C  
ATOM    126  O   ALA A  27      -5.809  54.154   4.602  1.00 72.15           O  
ATOM    127  CB  ALA A  27      -7.464  51.793   6.287  1.00 61.00           C  
ATOM    128  N   ALA A  28      -7.146  54.938  10.214  1.00 75.61           N  
ATOM    129  CA  ALA A  28      -5.786  54.427  10.004  1.00 75.21           C  
ATOM    130  C   ALA A  28      -5.518  53.781   8.645  1.00 73.11           C  
ATOM    131  O   ALA A  28      -4.360  53.601   8.258  1.00 75.35           O  
ATOM    132  N   ALA A  29      -9.995  54.900   9.697  1.00 80.60           N  
ATOM    133  CA  ALA A  29      -9.176  56.114   9.693  1.00 79.36           C  
ATOM    134  C   ALA A  29      -7.740  55.780   9.372  1.00 77.71           C  
ATOM    135  O   ALA A  29      -7.207  56.223   8.354  1.00 79.80           O  
ATOM    136  N   ALA A  30     -10.680  52.136   9.974  1.00 68.15           N  
ATOM    137  CA  ALA A  30     -11.040  52.912   8.797  1.00 70.17           C  
ATOM    138  C   ALA A  30     -10.146  54.136   8.619  1.00 77.31           C  
ATOM    139  O   ALA A  30      -9.592  54.365   7.534  1.00 81.98           O  
ATOM    140  CB  ALA A  30     -12.486  53.366   8.885  1.00 67.59           C  
ATOM    141  N   ALA A  31      -9.791  51.504  12.571  1.00 84.78           N  
ATOM    142  CA  ALA A  31      -9.188  50.876  11.395  1.00 82.80           C  
ATOM    143  C   ALA A  31      -9.444  51.691  10.145  1.00 76.45           C  
ATOM    144  O   ALA A  31      -8.530  51.913   9.355  1.00 76.25           O  
ATOM    145  CB  ALA A  31      -9.717  49.466  11.185  1.00 81.53           C  
ATOM    146  N   ALA A  32     -11.485  53.249  14.058  1.00 88.10           N  
ATOM    147  CA  ALA A  32     -10.033  53.239  14.263  1.00 88.48           C  
ATOM    148  C   ALA A  32      -9.260  52.593  13.117  1.00 86.25           C  
ATOM    149  O   ALA A  32      -8.188  53.067  12.755  1.00 90.90           O  
ATOM    150  CB  ALA A  32      -9.664  52.518  15.556  1.00 92.14           C  
ATOM    151  N   ALA A  33     -13.968  52.120  13.342  1.00 91.90           N  
ATOM    152  CA  ALA A  33     -13.558  53.413  12.773  1.00 91.90           C  
ATOM    153  C   ALA A  33     -12.043  53.633  12.909  1.00 89.11           C  
ATOM    154  O   ALA A  33     -11.398  54.100  11.970  1.00 83.38           O  
ATOM    155  CB  ALA A  33     -14.366  54.575  13.382  1.00 88.96           C  
ATOM    156  N   ALA A  34     -13.896  49.748  14.975  1.00 87.47           N  
ATOM    157  CA  ALA A  34     -14.087  49.721  13.520  1.00 86.78           C  
ATOM    158  C   ALA A  34     -13.504  50.959  12.879  1.00 89.36           C  
ATOM    159  O   ALA A  34     -12.679  50.860  11.969  1.00 90.48           O  
ATOM    160  CB  ALA A  34     -15.569  49.668  13.136  1.00 86.66           C  
ATOM    161  N   ALA A  35     -13.550  50.695  17.591  1.00 89.35           N  
ATOM    162  CA  ALA A  35     -12.513  49.840  17.010  1.00 90.77           C  
ATOM    163  C   ALA A  35     -12.710  49.495  15.525  1.00 89.57           C  
ATOM    164  O   ALA A  35     -11.787  49.001  14.888  1.00 90.71           O  
ATOM    165  CB  ALA A  35     -12.382  48.551  17.827  1.00 92.76           C  
ATOM    166  N   ALA A  36     -16.091  51.946  17.865  1.00 84.15           N  
ATOM    167  CA  ALA A  36     -14.839  52.725  17.945  1.00 83.59           C  
ATOM    168  C   ALA A  36     -13.671  51.977  17.276  1.00 85.46           C  
ATOM    169  O   ALA A  36     -12.921  52.545  16.480  1.00 82.35           O  
ATOM    170  CB  ALA A  36     -14.481  53.094  19.392  1.00 78.68           C  
TER     171      ALA A  36                                                      
ENDMDL                                                                          
END"""


class SearchModelTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.path.isdir(swamp.ENSEMBLE_DIR):
            os.makedirs(swamp.ENSEMBLE_DIR)
        if not os.path.isdir(swamp.IDEALHELICES_DIR):
            os.makedirs(swamp.IDEALHELICES_DIR)
        for fname in [os.path.join(swamp.ENSEMBLE_DIR, 'centroid_5.pdb'),
                      os.path.join(swamp.ENSEMBLE_DIR, 'ensemble_3.pdb'),
                      os.path.join(swamp.IDEALHELICES_DIR, 'ensemble_20_nativebfact_homogenous.pdb')]:
            with open(fname, 'w') as fhandle:
                fhandle.write(ENSEMBLE_SIDECHAINS)
            compress(fname)

    def test_1(self):

        searchmodel = SearchModel(id='test', ensemble_code='3', ermsd=0.6, nsearch=3, disable_check=False, mod='unmod',
                                  workdir=os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test'), model='ensemble')
        self.assertFalse(searchmodel.error)
        self.assertTrue(os.path.isdir(searchmodel.workdir))
        self.assertEqual(os.path.join(swamp.ENSEMBLE_DIR, 'ensemble_3.pdb.gz'), searchmodel.gzfile)
        self.addCleanup(remove, os.path.join(swamp.ENSEMBLE_DIR, 'ensemble_3.pdb.gz'))
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'ensemble_3.pdb'),
                         searchmodel.pdbfname)
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'ensemble_3.pdb'))
        self.assertTrue(os.path.isfile(searchmodel.pdbfname))
        self.assertDictEqual({'id': 'test', 'ermsd': 0.6, 'nsearch': 3, 'disable_check': False,
                              'pdbfile': os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'ensemble_3.pdb')},
                             searchmodel.phaser_info)

    @unittest.skipIf('THIS_IS_TRAVIS' in os.environ, "not implemented in Travis CI")
    def test_2(self):

        searchmodel = SearchModel(id='test', ensemble_code='5', ermsd=0.6, nsearch=3, disable_check=False,
                                  workdir=os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test'),
                                  mod='polyala', model='centroid')

        self.assertEqual(os.path.join(swamp.ENSEMBLE_DIR, 'centroid_5.pdb.gz'), searchmodel.gzfile)
        self.addCleanup(remove, os.path.join(swamp.ENSEMBLE_DIR, 'centroid_5.pdb.gz'))
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'centroid_5.pdb'))
        self.assertFalse(searchmodel.error)
        self.assertTrue(os.path.isfile(searchmodel.pdbfname))
        with open(searchmodel.pdbfname, 'r') as fhandle:
            contents = fhandle.readlines()
        self.assertListEqual([x.rstrip() for x in contents], [x.rstrip() for x in ENSEMBLE_POLYALA.split('\n')])

        self.assertListEqual([os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'models', 'model_1.pdb'),
                              os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'models', 'model_2.pdb')],
                             searchmodel.model_list)
        self.assertListEqual([os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'model_1_polyala.pdb'),
                              os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'model_2_polyala.pdb')],
                             searchmodel.modified_model_list)

        self.assertFalse(searchmodel.error)
        os.remove(searchmodel.modified_pdbfname)
        searchmodel._check_output()
        self.assertTrue(searchmodel.error)

    def test_3(self):
        searchmodel = SearchModel(id='test', ensemble_code='idealhelix', ermsd=0.6, nsearch=3, disable_check=False,
                                  workdir=os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test'),
                                  mod='unmod', model='centroid')

        self.assertEqual(os.path.join(swamp.IDEALHELICES_DIR, 'ensemble_20_nativebfact_homogenous.pdb.gz'),
                         searchmodel.idealhelix_fname)

        self.assertEqual(searchmodel.idealhelix_fname, searchmodel.gzfile)
        self.assertEqual(searchmodel.pdbfname, os.path.join(os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test',
                                                                         'idealhelix.pdb')))

    def test_4(self):

        fname_1 = create_tempfile(TEST_FULLSIZE_PDB_SIDECHAINS)
        self.addCleanup(remove, fname_1)
        fname_2 = get_tempfile()
        self.addCleanup(remove, fname_2)
        SearchModel.truncate_polyALA(fname_1, fname_2)
        self.assertTrue(os.path.isfile(fname_2))
        with open(fname_2, 'r') as fhandle:
            contents = fhandle.readlines()
        self.assertListEqual([x.rstrip() for x in TEST_FULLSIZE_PDB_POLYALA.split('\n')],
                             [x.rstrip() for x in contents])

    def test_5(self):
        fname_1 = create_tempfile(ENSEMBLE_POLYALA)
        self.addCleanup(remove, fname_1)
        fname_2 = create_tempfile(TEST_FULLSIZE_PDB_SIDECHAINS)
        self.addCleanup(remove, fname_2)
        SearchModel.transfer_flags_pdb(pdb_file=fname_1, pdb_ref=fname_2)
        with open(fname_1, 'r') as fhandle:
            contents = fhandle.readlines()
        self.assertListEqual([x.rstrip() for x in ENSEMBLE_POLYALA_WITH_FLAGS.split('\n')],
                             [x.rstrip() for x in contents])

    def test_6(self):
        fname_1 = create_tempfile(ENSEMBLE_POLYALA)
        self.addCleanup(remove, fname_1)
        os.mkdir(os.path.join(os.environ['CCP4_SCR'], 'models'))
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'models'))
        models = SearchModel.split_models(fname_1, directory=os.path.join(os.environ['CCP4_SCR'], 'models'),
                                          strip_hetatm=True)
        self.assertListEqual([os.path.join(os.environ['CCP4_SCR'], 'models', 'model_1.pdb'),
                              os.path.join(os.environ['CCP4_SCR'], 'models', 'model_2.pdb')],
                             models)
        self.assertListEqual(['model_1.pdb', 'model_2.pdb'],
                             sorted(os.listdir(os.path.join(os.environ['CCP4_SCR'], 'models'))))

    def test_7(self):
        searchmodel = SearchModel(id='test', ensemble_code='3', ermsd=0.6, nsearch=3, disable_check=False, mod='unmod',
                                  workdir=os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test'), model='ensemble')
        self.assertListEqual([os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'models', 'model_1.pdb'),
                              os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'models', 'model_2.pdb')],
                             searchmodel.model_list)

    def test_8(self):

        searchmodel = SearchModel(id='test', ensemble_code='999', ermsd=0.6, nsearch=3, disable_check=False,
                                  workdir=os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test'),
                                  mod='unmod', model='ensemble')

        self.assertTrue(searchmodel.error)
        self.addCleanup(remove, swamp.LIBRARY)
        self.addCleanup(remove, searchmodel.workdir)
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'models'), searchmodel.model_dir)
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'searchmodel_test', 'unmod_ensemble_999.pdb'),
                         searchmodel.modified_pdbfname)
