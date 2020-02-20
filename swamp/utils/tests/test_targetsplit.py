import os
import unittest
from swamp.utils import create_tempfile
from swamp.utils.targetsplit import TargetSplit


class TargetSplitTestCase(unittest.TestCase):

    def test_1(self):
        topcons_contents = """TOPCONS predicted topology:
iiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiii
"""
        pdb_contents = """CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2      8          
REMARK 465                                                                      
REMARK 465 MISSING RESIDUES                                                     
REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       
REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               
REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                
REMARK 465                                                                      
REMARK 465   M RES C SSSEQI                                                     
REMARK 465     MET A    -4                                                      
REMARK 465     VAL A    -3                                                      
REMARK 465     ALA A    -2                                                      
REMARK 465     ALA A    -1                                                      
REMARK 465     SER A     0                                                      
REMARK 465     MET A     1                                                      
REMARK 465     GLY A    98                                                      
REMARK 465     LYS A    99                                                      
REMARK 465     HIS A   212                                                      
REMARK 465     LYS A   215                                                      
ATOM    760  N   VAL A 100      17.668  61.385  96.142  1.00 36.12           N  
ANISOU  760  N   VAL A 100     4189   5832   3703    370    -20     96       N  
ATOM    761  CA  VAL A 100      16.510  62.175  95.720  1.00 34.76           C  
ANISOU  761  CA  VAL A 100     3981   5676   3550    300     62     84       C  
ATOM    762  C   VAL A 100      16.924  63.214  94.641  1.00 39.15           C  
ANISOU  762  C   VAL A 100     4461   6274   4139    307     77     -9       C  
ATOM    763  O   VAL A 100      16.205  63.379  93.656  1.00 38.11           O  
ANISOU  763  O   VAL A 100     4288   6134   4059    275    108    -15       O  
ATOM    764  CB  VAL A 100      15.715  62.769  96.916  1.00 37.75           C  
ANISOU  764  CB  VAL A 100     4379   6111   3852    257    129    130       C  
ATOM    765  CG1 VAL A 100      14.623  63.727  96.450  1.00 36.89           C  
ANISOU  765  CG1 VAL A 100     4216   6025   3776    215    217    110       C  
ATOM    766  CG2 VAL A 100      15.112  61.661  97.786  1.00 38.05           C  
ANISOU  766  CG2 VAL A 100     4485   6113   3858    228    124    244       C  
ATOM    767  N   GLY A 101      18.105  63.825  94.809  1.00 36.09           N  
ANISOU  767  N   GLY A 101     4052   5944   3718    343     50    -70       N  
ATOM    768  CA  GLY A 101      18.670  64.791  93.867  1.00 34.68           C  
ANISOU  768  CA  GLY A 101     3805   5805   3566    340     63   -145       C  
ATOM    769  C   GLY A 101      18.998  64.193  92.514  1.00 37.41           C  
ANISOU  769  C   GLY A 101     4110   6137   3967    361     26   -177       C  
ATOM    770  O   GLY A 101      18.818  64.843  91.481  1.00 35.74           O  
ANISOU  770  O   GLY A 101     3843   5954   3784    335     57   -198       O  
ATOM    771  N   VAL A 102      19.463  62.931  92.513  1.00 34.96           N  
ANISOU  771  N   VAL A 102     3830   5784   3671    410    -36   -177       N  
ATOM    772  CA  VAL A 102      19.819  62.187  91.297  1.00 34.18           C  
ANISOU  772  CA  VAL A 102     3699   5666   3623    436    -67   -233       C  
ATOM    773  C   VAL A 102      18.531  61.710  90.593  1.00 37.41           C  
ANISOU  773  C   VAL A 102     4118   6023   4073    373    -37   -212       C  
ATOM    774  O   VAL A 102      18.409  61.831  89.370  1.00 35.53           O  
ANISOU  774  O   VAL A 102     3822   5829   3850    347    -28   -263       O  
ATOM    775  CB  VAL A 102      20.820  61.047  91.624  1.00 38.45           C  
ANISOU  775  CB  VAL A 102     4268   6161   4180    528   -137   -249       C  
ATOM    776  CG1 VAL A 102      21.126  60.185  90.399  1.00 38.40           C  
ANISOU  776  CG1 VAL A 102     4237   6117   4236    561   -155   -331       C  
ATOM    777  CG2 VAL A 102      22.111  61.608  92.229  1.00 37.89           C  
ANISOU  777  CG2 VAL A 102     4155   6188   4054    582   -172   -273       C  
ATOM    778  N   ILE A 103      17.542  61.236  91.381  1.00 34.36           N  
ANISOU  778  N   ILE A 103     3794   5569   3692    337    -19   -133       N  
ATOM    779  CA  ILE A 103      16.260  60.794  90.844  1.00 33.66           C  
ANISOU  779  CA  ILE A 103     3704   5449   3636    259     10   -110       C  
ATOM    780  C   ILE A 103      15.544  61.966  90.187  1.00 37.70           C  
ANISOU  780  C   ILE A 103     4134   6061   4131    214     61   -102       C  
ATOM    781  O   ILE A 103      15.031  61.813  89.070  1.00 37.84           O  
ANISOU  781  O   ILE A 103     4097   6118   4163    169     63   -129       O  
ATOM    782  CB  ILE A 103      15.417  60.020  91.896  1.00 37.07           C  
ANISOU  782  CB  ILE A 103     4214   5798   4074    220     22    -17       C  
ATOM    783  CG1 ILE A 103      16.062  58.633  92.170  1.00 37.34           C  
ANISOU  783  CG1 ILE A 103     4330   5703   4156    266    -32    -15       C  
ATOM    784  CG2 ILE A 103      13.920  59.876  91.451  1.00 37.66           C  
ANISOU  784  CG2 ILE A 103     4258   5883   4167    115     67     14       C  
ATOM    785  CD1 ILE A 103      15.598  57.949  93.432  1.00 47.42           C  
ANISOU  785  CD1 ILE A 103     5694   6900   5425    246    -28    105       C  
ATOM    786  N   LEU A 104      15.594  63.153  90.831  1.00 33.68           N  
ANISOU  786  N   LEU A 104     3611   5596   3591    231    101    -71       N  
ATOM    787  CA  LEU A 104      14.977  64.376  90.307  1.00 33.36           C  
ANISOU  787  CA  LEU A 104     3499   5623   3552    211    157    -47       C  
ATOM    788  C   LEU A 104      15.511  64.746  88.917  1.00 34.07           C  
ANISOU  788  C   LEU A 104     3518   5780   3648    214    139    -90       C  
ATOM    789  O   LEU A 104      14.708  65.012  88.027  1.00 32.90           O  
ANISOU  789  O   LEU A 104     3303   5692   3505    182    158    -56       O  
ATOM    790  CB  LEU A 104      15.136  65.530  91.310  1.00 33.69           C  
ANISOU  790  CB  LEU A 104     3558   5667   3576    234    207    -34       C  
ATOM    791  CG  LEU A 104      14.360  66.816  91.054  1.00 38.02           C  
ANISOU  791  CG  LEU A 104     4051   6243   4152    231    283      6       C  
ATOM    792  CD1 LEU A 104      12.849  66.546  90.864  1.00 37.78           C  
ANISOU  792  CD1 LEU A 104     3981   6236   4138    203    317     79       C  
ATOM    793  CD2 LEU A 104      14.564  67.790  92.218  1.00 40.00           C  
ANISOU  793  CD2 LEU A 104     4341   6468   4390    248    340    -14       C  
ATOM    794  N   VAL A 105      16.858  64.715  88.727  1.00 31.14           N  
ANISOU  794  N   VAL A 105     3147   5419   3264    251    101   -159       N  
ATOM    795  CA  VAL A 105      17.526  64.969  87.443  1.00 30.96           C  
ANISOU  795  CA  VAL A 105     3053   5478   3231    250     87   -207       C  
ATOM    796  C   VAL A 105      17.042  63.944  86.411  1.00 35.03           C  
ANISOU  796  C   VAL A 105     3543   6022   3746    218     59   -243       C  
ATOM    797  O   VAL A 105      16.709  64.332  85.295  1.00 35.34           O  
ANISOU  797  O   VAL A 105     3507   6161   3761    183     71   -233       O  
ATOM    798  CB  VAL A 105      19.074  64.920  87.574  1.00 35.01           C  
ANISOU  798  CB  VAL A 105     3565   6008   3728    296     51   -284       C  
ATOM    799  CG1 VAL A 105      19.755  64.761  86.208  1.00 34.91           C  
ANISOU  799  CG1 VAL A 105     3479   6090   3696    293     33   -351       C  
ATOM    800  CG2 VAL A 105      19.596  66.143  88.285  1.00 34.55           C  
ANISOU  800  CG2 VAL A 105     3509   5958   3662    295     82   -268       C  
ATOM    801  N   GLY A 106      17.024  62.660  86.802  1.00 31.90           N  
ANISOU  801  N   GLY A 106     3210   5538   3374    226     26   -284       N  
ATOM    802  CA  GLY A 106      16.603  61.553  85.954  1.00 32.79           C  
ANISOU  802  CA  GLY A 106     3316   5644   3498    184      4   -348       C  
ATOM    803  C   GLY A 106      15.165  61.656  85.488  1.00 38.66           C  
ANISOU  803  C   GLY A 106     4015   6444   4230     96     28   -295       C  
ATOM    804  O   GLY A 106      14.840  61.246  84.373  1.00 39.30           O  
ANISOU  804  O   GLY A 106     4041   6603   4287     40     16   -354       O  
ATOM    805  N   CYS A 107      14.292  62.202  86.336  1.00 34.66           N  
ANISOU  805  N   CYS A 107     3520   5918   3732     83     64   -190       N  
ATOM    806  CA  CYS A 107      12.871  62.327  86.029  1.00 33.96           C  
ANISOU  806  CA  CYS A 107     3372   5898   3635     11     89   -126       C  
ATOM    807  C   CYS A 107      12.559  63.546  85.180  1.00 36.79           C  
ANISOU  807  C   CYS A 107     3622   6399   3958     17    112    -65       C  
ATOM    808  O   CYS A 107      11.462  63.650  84.641  1.00 34.93           O  
ANISOU  808  O   CYS A 107     3306   6262   3703    -34    121    -13       O  
ATOM    809  CB  CYS A 107      12.047  62.300  87.309  1.00 34.50           C  
ANISOU  809  CB  CYS A 107     3487   5894   3726      0    126    -44       C  
ATOM    810  SG  CYS A 107      12.085  60.707  88.159  1.00 39.40           S  
ANISOU  810  SG  CYS A 107     4227   6357   4386    -36     99    -71       S  
ATOM    811  N   CYS A 108      13.515  64.471  85.058  1.00 35.29           N  
ANISOU  811  N   CYS A 108     3424   6224   3761     77    122    -62       N  
ATOM    812  CA  CYS A 108      13.303  65.682  84.256  1.00 35.81           C  
ANISOU  812  CA  CYS A 108     3398   6402   3805     88    148     20       C  
ATOM    813  C   CYS A 108      13.248  65.386  82.748  1.00 39.76           C  
ANISOU  813  C   CYS A 108     3808   7060   4239     38    113     -8       C  
ATOM    814  O   CYS A 108      13.805  64.369  82.295  1.00 39.23           O  
ANISOU  814  O   CYS A 108     3760   7001   4146      8     73   -131       O  
ATOM    815  CB  CYS A 108      14.373  66.725  84.577  1.00 35.44           C  
ANISOU  815  CB  CYS A 108     3377   6313   3777    143    174     28       C  
ATOM    816  SG  CYS A 108      14.063  67.645  86.106  1.00 38.80           S  
ANISOU  816  SG  CYS A 108     3867   6613   4263    189    240     86       S  
ATOM    817  N   PRO A 109      12.626  66.278  81.941  1.00 36.74           N  
ANISOU  817  N   PRO A 109     3324   6811   3825     33    128    102       N  
ATOM    818  CA  PRO A 109      12.651  66.072  80.487  1.00 36.54           C  
ANISOU  818  CA  PRO A 109     3203   6973   3707    -19     92     82       C  
ATOM    819  C   PRO A 109      14.051  66.355  79.917  1.00 39.90           C  
ANISOU  819  C   PRO A 109     3632   7433   4096     -2     88     25       C  
ATOM    820  O   PRO A 109      14.984  66.720  80.652  1.00 38.31           O  
ANISOU  820  O   PRO A 109     3500   7110   3946     47    109      3       O  
ATOM    821  CB  PRO A 109      11.626  67.086  79.984  1.00 38.91           C  
ANISOU  821  CB  PRO A 109     3395   7399   3989     -7    112    256       C  
ATOM    822  CG  PRO A 109      11.709  68.211  80.953  1.00 43.16           C  
ANISOU  822  CG  PRO A 109     3985   7792   4623     78    173    354       C  
ATOM    823  CD  PRO A 109      11.921  67.533  82.291  1.00 38.62           C  
ANISOU  823  CD  PRO A 109     3527   7039   4109     85    182    254       C  
ATOM    824  N   GLY A 110      14.178  66.209  78.606  1.00 36.96           N  
ANISOU  824  N   GLY A 110     3172   7251   3622    -51     63      1       N  
ATOM    825  CA  GLY A 110      15.412  66.500  77.896  1.00 36.79           C  
ANISOU  825  CA  GLY A 110     3125   7311   3543    -47     66    -43       C  
ATOM    826  C   GLY A 110      15.754  67.977  77.871  1.00 39.68           C  
ANISOU  826  C   GLY A 110     3468   7678   3930    -11    108    116       C  
ATOM    827  O   GLY A 110      14.932  68.832  78.237  1.00 37.90           O  
ANISOU  827  O   GLY A 110     3235   7404   3761     20    137    269       O  
ATOM    828  N   GLY A 111      16.997  68.255  77.492  1.00 37.57           N  
ANISOU  828  N   GLY A 111     3192   7452   3630    -14    120     73       N  
ATOM    829  CA  GLY A 111      17.524  69.610  77.421  1.00 38.32           C  
ANISOU  829  CA  GLY A 111     3273   7535   3751     -3    166    207       C  
ATOM    830  C   GLY A 111      17.587  70.139  76.005  1.00 45.38           C  
ANISOU  830  C   GLY A 111     4057   8656   4530    -51    169    312       C  
ATOM    831  O   GLY A 111      17.838  69.373  75.069  1.00 44.83           O  
ANISOU  831  O   GLY A 111     3925   8771   4339    -97    138    211       O  
ATOM    832  N   THR A 112      17.360  71.460  75.849  1.00 44.03           N  
ANISOU  832  N   THR A 112     3865   8469   4397    -38    211    517       N  
ATOM    833  CA  THR A 112      17.408  72.178  74.572  1.00 46.06           C  
ANISOU  833  CA  THR A 112     4021   8928   4552    -78    221    676       C  
ATOM    834  C   THR A 112      18.755  71.978  73.857  1.00 50.30           C  
ANISOU  834  C   THR A 112     4519   9608   4985   -142    226    577       C  
ATOM    835  O   THR A 112      18.772  71.719  72.654  1.00 51.95           O  
ANISOU  835  O   THR A 112     4629  10072   5039   -195    206    587       O  
ATOM    836  CB  THR A 112      17.000  73.648  74.806  1.00 62.92           C  
ANISOU  836  CB  THR A 112     6172  10938   6798    -35    277    915       C  
ATOM    837  OG1 THR A 112      15.573  73.736  74.755  1.00 69.65           O  
ANISOU  837  OG1 THR A 112     6981  11815   7666     17    260   1043       O  
ATOM    838  CG2 THR A 112      17.610  74.618  73.802  1.00 64.44           C  
ANISOU  838  CG2 THR A 112     6303  11251   6931    -82    311   1083       C  
ATOM    839  N   ALA A 113      19.873  72.042  74.608  1.00 44.39           N  
ANISOU  839  N   ALA A 113     3836   8720   4311   -139    253    467       N  
ATOM    840  CA  ALA A 113      21.222  71.889  74.071  1.00 44.43           C  
ANISOU  840  CA  ALA A 113     3795   8854   4234   -190    266    364       C  
ATOM    841  C   ALA A 113      21.413  70.625  73.231  1.00 47.29           C  
ANISOU  841  C   ALA A 113     4088   9428   4452   -211    225    186       C  
ATOM    842  O   ALA A 113      22.263  70.622  72.346  1.00 49.14           O  
ANISOU  842  O   ALA A 113     4242   9859   4571   -263    243    147       O  
ATOM    843  CB  ALA A 113      22.253  71.957  75.192  1.00 44.42           C  
ANISOU  843  CB  ALA A 113     3868   8672   4339   -172    285    250       C  
ATOM    844  N   SER A 114      20.592  69.580  73.459  1.00 42.02           N  
ANISOU  844  N   SER A 114     3449   8728   3789   -181    179     77       N  
ATOM    845  CA  SER A 114      20.643  68.338  72.676  1.00 42.09           C  
ANISOU  845  CA  SER A 114     3406   8907   3678   -208    146   -113       C  
ATOM    846  C   SER A 114      20.318  68.595  71.190  1.00 48.23           C  
ANISOU  846  C   SER A 114     4059   9997   4271   -285    144    -27       C  
ATOM    847  O   SER A 114      20.812  67.864  70.334  1.00 49.03           O  
ANISOU  847  O   SER A 114     4094  10295   4241   -326    141   -189       O  
ATOM    848  CB  SER A 114      19.694  67.294  73.251  1.00 42.07           C  
ANISOU  848  CB  SER A 114     3468   8782   3735   -182    104   -217       C  
ATOM    849  OG  SER A 114      18.338  67.653  73.046  1.00 40.48           O  
ANISOU  849  OG  SER A 114     3237   8626   3517   -204     85    -61       O  
ATOM    850  N   ASN A 115      19.509  69.644  70.892  1.00 45.56           N  
ANISOU  850  N   ASN A 115     3683   9708   3919   -299    148    228       N  
ATOM    851  CA  ASN A 115      19.145  70.032  69.522  1.00 47.16           C  
ANISOU  851  CA  ASN A 115     3760  10222   3938   -367    140    364       C  
ATOM    852  C   ASN A 115      20.388  70.481  68.749  1.00 52.79           C  
ANISOU  852  C   ASN A 115     4408  11108   4543   -421    186    376       C  
ATOM    853  O   ASN A 115      20.575  70.064  67.608  1.00 53.35           O  
ANISOU  853  O   ASN A 115     4376  11476   4418   -488    179    307       O  
ATOM    854  CB  ASN A 115      18.077  71.140  69.504  1.00 44.25           C  
ANISOU  854  CB  ASN A 115     3369   9833   3610   -341    138    667       C  
ATOM    855  CG  ASN A 115      16.800  70.842  70.247  1.00 49.80           C  
ANISOU  855  CG  ASN A 115     4113  10392   4415   -288    102    684       C  
ATOM    856  OD1 ASN A 115      16.532  69.712  70.684  1.00 42.29           O  
ANISOU  856  OD1 ASN A 115     3203   9378   3486   -290     71    477       O  
ATOM    857  ND2 ASN A 115      15.969  71.870  70.401  1.00 36.54           N  
ANISOU  857  ND2 ASN A 115     2419   8656   2807   -236    113    942       N  
ATOM    858  N   VAL A 116      21.239  71.306  69.392  1.00 50.44           N  
ANISOU  858  N   VAL A 116     4165  10635   4366   -402    236    447       N  
ATOM    859  CA  VAL A 116      22.505  71.841  68.853  1.00 52.16           C  
ANISOU  859  CA  VAL A 116     4326  10981   4511   -463    290    468       C  
ATOM    860  C   VAL A 116      23.544  70.715  68.752  1.00 54.64           C  
ANISOU  860  C   VAL A 116     4613  11384   4763   -465    293    166       C  
ATOM    861  O   VAL A 116      24.294  70.673  67.783  1.00 55.35           O  
ANISOU  861  O   VAL A 116     4601  11734   4694   -529    323    124       O  
ATOM    862  CB  VAL A 116      23.071  73.022  69.706  1.00 56.89           C  
ANISOU  862  CB  VAL A 116     4998  11337   5279   -457    344    609       C  
ATOM    863  CG1 VAL A 116      24.071  73.845  68.896  1.00 58.16           C  
ANISOU  863  CG1 VAL A 116     5080  11672   5348   -553    403    725       C  
ATOM    864  CG2 VAL A 116      21.955  73.923  70.239  1.00 56.78           C  
ANISOU  864  CG2 VAL A 116     5052  11120   5402   -409    344    840       C  
ATOM    865  N   MET A 117      23.610  69.824  69.766  1.00 49.13           N  
ANISOU  865  N   MET A 117     4004  10473   4191   -388    266    -33       N  
ATOM    866  CA  MET A 117      24.552  68.701  69.767  1.00 48.45           C  
ANISOU  866  CA  MET A 117     3899  10430   4078   -358    268   -314       C  
ATOM    867  C   MET A 117      24.254  67.711  68.648  1.00 52.17           C  
ANISOU  867  C   MET A 117     4293  11152   4378   -393    252   -474       C  
ATOM    868  O   MET A 117      25.194  67.215  68.028  1.00 52.50           O  
ANISOU  868  O   MET A 117     4257  11371   4319   -405    284   -646       O  
ATOM    869  CB  MET A 117      24.629  68.000  71.128  1.00 49.33           C  
ANISOU  869  CB  MET A 117     4128  10247   4370   -260    239   -453       C  
ATOM    870  CG  MET A 117      25.150  68.879  72.238  1.00 52.66           C  
ANISOU  870  CG  MET A 117     4613  10460   4935   -235    258   -352       C  
ATOM    871  SD  MET A 117      26.786  69.586  71.949  1.00 59.01           S  
ANISOU  871  SD  MET A 117     5326  11401   5696   -287    317   -356       S  
ATOM    872  CE  MET A 117      26.451  71.275  72.409  1.00 55.97           C  
ANISOU  872  CE  MET A 117     4993  10867   5406   -350    352    -71       C  
ATOM    873  N   THR A 118      22.954  67.453  68.369  1.00 48.22           N  
ANISOU  873  N   THR A 118     3801  10684   3838   -416    208   -425       N  
ATOM    874  CA  THR A 118      22.510  66.589  67.261  1.00 49.35           C  
ANISOU  874  CA  THR A 118     3864  11087   3801   -477    189   -573       C  
ATOM    875  C   THR A 118      22.875  67.211  65.915  1.00 54.31           C  
ANISOU  875  C   THR A 118     4351  12081   4202   -570    221   -475       C  
ATOM    876  O   THR A 118      23.249  66.476  64.997  1.00 55.68           O  
ANISOU  876  O   THR A 118     4443  12501   4211   -616    237   -677       O  
ATOM    877  CB  THR A 118      21.012  66.325  67.336  1.00 50.39           C  
ANISOU  877  CB  THR A 118     4022  11182   3942   -497    131   -513       C  
ATOM    878  OG1 THR A 118      20.724  65.778  68.612  1.00 43.18           O  
ANISOU  878  OG1 THR A 118     3238   9934   3233   -419    110   -590       O  
ATOM    879  CG2 THR A 118      20.537  65.369  66.270  1.00 47.96           C  
ANISOU  879  CG2 THR A 118     3636  11134   3454   -578    107   -697       C  
ATOM    880  N   TYR A 119      22.760  68.554  65.795  1.00 51.06           N  
ANISOU  880  N   TYR A 119     3914  11704   3783   -598    236   -166       N  
ATOM    881  CA  TYR A 119      23.117  69.271  64.566  1.00 53.65           C  
ANISOU  881  CA  TYR A 119     4114  12372   3900   -691    271    -16       C  
ATOM    882  C   TYR A 119      24.628  69.129  64.296  1.00 59.16           C  
ANISOU  882  C   TYR A 119     4757  13180   4541   -710    338   -173       C  
ATOM    883  O   TYR A 119      25.034  68.875  63.162  1.00 59.37           O  
ANISOU  883  O   TYR A 119     4665  13547   4347   -783    366   -255       O  
ATOM    884  CB  TYR A 119      22.679  70.751  64.632  1.00 55.51           C  
ANISOU  884  CB  TYR A 119     4355  12551   4185   -703    279    366       C  
ATOM    885  CG  TYR A 119      23.166  71.571  63.454  1.00 60.61           C  
ANISOU  885  CG  TYR A 119     4880  13516   4633   -800    322    559       C  
ATOM    886  CD1 TYR A 119      22.548  71.478  62.209  1.00 63.65           C  
ANISOU  886  CD1 TYR A 119     5144  14271   4768   -873    293    637       C  
ATOM    887  CD2 TYR A 119      24.279  72.401  63.569  1.00 62.92           C  
ANISOU  887  CD2 TYR A 119     5170  13763   4972   -833    393    656       C  
ATOM    888  CE1 TYR A 119      23.026  72.193  61.109  1.00 66.91           C  
ANISOU  888  CE1 TYR A 119     5441  15003   4978   -968    334    824       C  
ATOM    889  CE2 TYR A 119      24.753  73.132  62.484  1.00 65.74           C  
ANISOU  889  CE2 TYR A 119     5416  14418   5143   -936    440    843       C  
ATOM    890  CZ  TYR A 119      24.127  73.022  61.255  1.00 76.19           C  
ANISOU  890  CZ  TYR A 119     6624  16110   6213   -999    411    933       C  
ATOM    891  OH  TYR A 119      24.601  73.748  60.192  1.00 84.36           O  
ANISOU  891  OH  TYR A 119     7548  17454   7051  -1105    458   1138       O  
ATOM    892  N   LEU A 120      25.443  69.243  65.355  1.00 56.13           N  
ANISOU  892  N   LEU A 120     4451  12528   4349   -644    363   -228       N  
ATOM    893  CA  LEU A 120      26.885  69.097  65.249  1.00 57.63           C  
ANISOU  893  CA  LEU A 120     4579  12807   4510   -648    422   -379       C  
ATOM    894  C   LEU A 120      27.279  67.646  64.978  1.00 60.25           C  
ANISOU  894  C   LEU A 120     4879  13229   4783   -599    423   -735       C  
ATOM    895  O   LEU A 120      28.129  67.411  64.125  1.00 61.78           O  
ANISOU  895  O   LEU A 120     4956  13698   4819   -639    477   -862       O  
ATOM    896  CB  LEU A 120      27.592  69.691  66.475  1.00 57.46           C  
ANISOU  896  CB  LEU A 120     4636  12495   4700   -601    440   -323       C  
ATOM    897  CG  LEU A 120      27.363  71.209  66.696  1.00 63.48           C  
ANISOU  897  CG  LEU A 120     5431  13154   5533   -661    461      9       C  
ATOM    898  CD1 LEU A 120      27.684  71.611  68.123  1.00 62.41           C  
ANISOU  898  CD1 LEU A 120     5407  12675   5630   -606    459     21       C  
ATOM    899  CD2 LEU A 120      28.144  72.063  65.682  1.00 68.31           C  
ANISOU  899  CD2 LEU A 120     5923  14043   5990   -782    529    157       C  
ATOM   1976  N   ALA A 270      17.395  71.560  89.366  1.00 35.38           N  
ANISOU 1976  N   ALA A 270     3622   5894   3926    228    397    -64       N  
ATOM   1977  CA  ALA A 270      17.973  72.447  90.379  1.00 35.08           C  
ANISOU 1977  CA  ALA A 270     3634   5791   3903    205    443   -136       C  
ATOM   1978  C   ALA A 270      16.885  73.370  90.924  1.00 38.15           C  
ANISOU 1978  C   ALA A 270     4048   6088   4360    230    542   -100       C  
ATOM   1979  O   ALA A 270      16.822  73.575  92.129  1.00 38.62           O  
ANISOU 1979  O   ALA A 270     4161   6109   4402    225    578   -172       O  
ATOM   1980  CB  ALA A 270      19.128  73.259  89.793  1.00 35.65           C  
ANISOU 1980  CB  ALA A 270     3680   5870   3995    157    449   -165       C  
ATOM   1981  N   ALA A 271      15.997  73.864  90.048  1.00 34.42           N  
ANISOU 1981  N   ALA A 271     3528   5596   3955    265    584     11       N  
ATOM   1982  CA  ALA A 271      14.866  74.726  90.424  1.00 34.94           C  
ANISOU 1982  CA  ALA A 271     3597   5575   4102    317    682     61       C  
ATOM   1983  C   ALA A 271      13.844  73.999  91.306  1.00 38.84           C  
ANISOU 1983  C   ALA A 271     4101   6101   4557    346    691     53       C  
ATOM   1984  O   ALA A 271      13.371  74.567  92.293  1.00 38.92           O  
ANISOU 1984  O   ALA A 271     4146   6047   4593    369    773      7       O  
ATOM   1985  CB  ALA A 271      14.180  75.270  89.179  1.00 35.46           C  
ANISOU 1985  CB  ALA A 271     3588   5648   4238    361    705    211       C  
ATOM   1986  N   LEU A 272      13.492  72.752  90.934  1.00 35.78           N  
ANISOU 1986  N   LEU A 272     3680   5809   4104    337    616     91       N  
ATOM   1987  CA  LEU A 272      12.516  71.931  91.663  1.00 35.50           C  
ANISOU 1987  CA  LEU A 272     3648   5811   4029    343    620    101       C  
ATOM   1988  C   LEU A 272      13.015  71.552  93.053  1.00 39.08           C  
ANISOU 1988  C   LEU A 272     4184   6252   4414    316    617      5       C  
ATOM   1989  O   LEU A 272      12.231  71.522  93.998  1.00 40.15           O  
ANISOU 1989  O   LEU A 272     4334   6392   4529    325    673      0       O  
ATOM   1990  CB  LEU A 272      12.119  70.681  90.845  1.00 34.66           C  
ANISOU 1990  CB  LEU A 272     3494   5793   3881    317    540    155       C  
ATOM   1991  CG  LEU A 272      11.233  70.919  89.606  1.00 38.06           C  
ANISOU 1991  CG  LEU A 272     3821   6290   4350    336    544    263       C  
ATOM   1992  CD1 LEU A 272      10.998  69.623  88.830  1.00 36.79           C  
ANISOU 1992  CD1 LEU A 272     3620   6223   4134    283    461    273       C  
ATOM   1993  CD2 LEU A 272       9.888  71.527  89.980  1.00 39.50           C  
ANISOU 1993  CD2 LEU A 272     3951   6476   4583    389    629    334       C  
ATOM   1994  N   ALA A 273      14.325  71.313  93.181  1.00 35.62           N  
ANISOU 1994  N   ALA A 273     3787   5818   3930    284    555    -67       N  
ATOM   1995  CA  ALA A 273      14.977  70.980  94.446  1.00 35.61           C  
ANISOU 1995  CA  ALA A 273     3851   5833   3846    260    535   -149       C  
ATOM   1996  C   ALA A 273      15.034  72.213  95.357  1.00 41.32           C  
ANISOU 1996  C   ALA A 273     4608   6511   4582    252    627   -233       C  
ATOM   1997  O   ALA A 273      14.832  72.086  96.566  1.00 41.96           O  
ANISOU 1997  O   ALA A 273     4729   6623   4590    240    654   -281       O  
ATOM   1998  CB  ALA A 273      16.385  70.473  94.176  1.00 35.68           C  
ANISOU 1998  CB  ALA A 273     3867   5876   3813    240    441   -195       C  
ATOM   1999  N   ALA A 274      15.306  73.400  94.773  1.00 37.99           N  
ANISOU 1999  N   ALA A 274     4172   6014   4250    253    679   -252       N  
ATOM   2000  CA  ALA A 274      15.361  74.677  95.495  1.00 38.62           C  
ANISOU 2000  CA  ALA A 274     4289   6011   4372    241    781   -348       C  
ATOM   2001  C   ALA A 274      13.976  75.054  96.014  1.00 40.49           C  
ANISOU 2001  C   ALA A 274     4522   6213   4650    298    887   -328       C  
ATOM   2002  O   ALA A 274      13.870  75.527  97.136  1.00 42.06           O  
ANISOU 2002  O   ALA A 274     4764   6399   4817    285    958   -435       O  
ATOM   2003  CB  ALA A 274      15.902  75.778  94.584  1.00 39.90           C  
ANISOU 2003  CB  ALA A 274     4439   6076   4645    227    814   -342       C  
ATOM   2004  N   ALA A 275      12.917  74.808  95.219  1.00 35.17           N  
ANISOU 2004  N   ALA A 275     3785   5545   4034    358    895   -198       N  
ATOM   2005  CA  ALA A 275      11.534  75.127  95.578  1.00 35.50           C  
ANISOU 2005  CA  ALA A 275     3792   5576   4119    425    993   -162       C  
ATOM   2006  C   ALA A 275      10.891  74.133  96.544  1.00 41.37           C  
ANISOU 2006  C   ALA A 275     4539   6428   4750    407    987   -168       C  
ATOM   2007  O   ALA A 275      10.161  74.562  97.441  1.00 41.44           O  
ANISOU 2007  O   ALA A 275     4550   6440   4754    435   1088   -219       O  
ATOM   2008  CB  ALA A 275      10.682  75.243  94.323  1.00 35.94           C  
ANISOU 2008  CB  ALA A 275     3756   5631   4267    491    994    -11       C  
ATOM   2009  N   HIS A 276      11.151  72.806  96.365  1.00 37.48           N  
ANISOU 2009  N   HIS A 276     4048   6022   4172    360    877   -116       N  
ATOM   2010  CA  HIS A 276      10.514  71.766  97.173  1.00 36.82           C  
ANISOU 2010  CA  HIS A 276     3970   6029   3991    331    867    -88       C  
ATOM   2011  C   HIS A 276      11.327  71.079  98.243  1.00 39.86           C  
ANISOU 2011  C   HIS A 276     4431   6468   4246    274    814   -146       C  
ATOM   2012  O   HIS A 276      10.741  70.470  99.133  1.00 40.60           O  
ANISOU 2012  O   HIS A 276     4537   6634   4256    250    835   -122       O  
ATOM   2013  CB  HIS A 276       9.784  70.774  96.275  1.00 37.80           C  
ANISOU 2013  CB  HIS A 276     4031   6199   4131    322    810     33       C  
ATOM   2014  CG  HIS A 276       8.828  71.460  95.363  1.00 41.63           C  
ANISOU 2014  CG  HIS A 276     4421   6676   4720    382    864    104       C  
ATOM   2015  ND1 HIS A 276       7.594  71.905  95.817  1.00 44.10           N  
ANISOU 2015  ND1 HIS A 276     4674   7023   5060    428    967    129       N  
ATOM   2016  CD2 HIS A 276       8.988  71.835  94.071  1.00 43.54           C  
ANISOU 2016  CD2 HIS A 276     4613   6894   5037    410    830    160       C  
ATOM   2017  CE1 HIS A 276       7.034  72.512  94.787  1.00 44.28           C  
ANISOU 2017  CE1 HIS A 276     4609   7036   5180    492    987    207       C  
ATOM   2018  NE2 HIS A 276       7.819  72.471  93.702  1.00 44.20           N  
ANISOU 2018  NE2 HIS A 276     4602   6997   5195    478    903    236       N   
"""
        conpred_contents = """PFRMAT RR
TARGET 536987
AUTHOR RaptorX-Contact
METHOD deep dilated residual networks (one variant of deep CNN). Consult jinboxu@gmail.com for details.
MODEL 1
MVAASMNILSKISSFIGKTFSLWAALFAAAAFFAPDTFKWAGPYIPWLLG
IIMFGMGLTLKPSDFDILFKHPKVVIIGVIAQFAIMPATAWCLSKLLNLP
AEIAVGVILVGCCPGGTASNVMTYLARGNVALSVAVTSVSTLTSPLLTPA
IFLMLAGEMLEIQAAGMLMSIVKMVLLPIVLGLIVHKVLGSKTEKLTDAL
PLVSVAAIVLIIGAVVGASKGKIMESGLLIFAVVVLHNGIGYLLGFFAAK
WTGLPYDAQKALTIEVGMQNSGLAAALAAAHFAAAPVVAVPGALFSVWHN
ISGSLLATYWAAKAGKHKKPLDRAGSENLYFQ
53 178 0 8 0.9999614
57 182 0 8 0.9999346
58 182 0 8 0.9999014
54 181 0 8 0.9998163
54 182 0 8 0.9997769
54 178 0 8 0.9996910
249 259 0 8 0.9989253
58 185 0 8 0.9979285
58 186 0 8 0.9977884
249 262 0 8 0.9974785
94 104 0 8 0.9972718
123 133 0 8 0.9972159
57 179 0 8 0.9963613
246 263 0 8 0.9962631
50 178 0 8 0.9946589
106 288 0 8 0.9932054
57 183 0 8 0.9925978
123 261 0 8 0.9922032
102 288 0 8 0.9917381
27 212 0 8 0.9908113
103 291 0 8 0.9907801
75 136 0 8 0.9905434
31 216 0 8 0.9904293
89 240 0 8 0.9902470
27 213 0 8 0.9900678
110 292 0 8 0.9887912
85 244 0 8 0.9886514
90 108 0 8 0.9883336
109 278 0 8 0.9877242
94 107 0 8 0.9875522
78 262 0 8 0.9875078
48 207 0 8 0.9874308
74 262 0 8 0.9874212
28 216 0 8 0.9870313
245 263 0 8 0.9866461
78 136 0 8 0.9865698
106 291 0 8 0.9861109
79 139 0 8 0.9859405
133 265 0 8 0.9857825
77 252 0 8 0.9857346
109 274 0 8 0.9857225
110 295 0 8 0.9855377
81 248 0 8 0.9851450
81 266 0 8 0.9848748
74 258 0 8 0.9841593
106 292 0 8 0.9837796
31 213 0 8 0.9835263
68 135 0 8 0.9834397
48 211 0 8 0.9833449
113 274 0 8 0.9828007
52 207 0 8 0.9818235
128 261 0 8 0.9814836
90 107 0 8 0.9814461
119 265 0 8 0.9814367
105 288 0 8 0.9791791
271 296 0 8 0.9788657
90 111 0 8 0.9781752
31 217 0 8 0.9776807
53 175 0 8 0.9772123
77 262 0 8 0.9764582
129 258 0 8 0.9764170
234 298 0 8 0.9763948
133 261 0 8 0.9759184
79 140 0 8 0.9759070
55 182 0 8 0.9758528
246 259 0 8 0.9756561
27 209 0 8 0.9746038
234 295 0 8 0.9741930
112 148 0 8 0.9737659
102 287 0 8 0.9732612
132 258 0 8 0.9728087
82 266 0 8 0.9718467
242 263 0 8 0.9710815
245 266 0 8 0.9700539
91 108 0 8 0.9698529
75 139 0 8 0.9698042
48 210 0 8 0.9697683
24 212 0 8 0.9695854
107 233 0 8 0.9683198
136 262 0 8 0.9669924
107 291 0 8 0.9663849
79 136 0 8 0.9657449
94 108 0 8 0.9650769
125 307 0 8 0.9650706
77 248 0 8 0.9650462
120 133 0 8 0.9647374
93 233 0 8 0.9635152
51 207 0 8 0.9634590
"""
        topcons_fname = create_tempfile(content=topcons_contents)
        self.addCleanup(os.remove, topcons_fname)
        conpred_fname = create_tempfile(content=conpred_contents)
        self.addCleanup(os.remove, conpred_fname)
        pdb_fname = create_tempfile(content=pdb_contents)
        self.addCleanup(os.remove, pdb_fname)
        targetsplit = TargetSplit(conpred=conpred_fname, sspred=topcons_fname, workdir=os.environ['CCP4_SCR'],
                                  pdb_benchmark=pdb_fname)
        targetsplit.split()

        self.assertEqual(16, len(targetsplit.subtargets_pdb['1_10'][0][0]))
        self.assertListEqual([15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 33, 34],
                             [x.seqid.num for x in targetsplit.subtargets_pdb['1_2'][0][0]])
        self.assertListEqual(
            [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (2, 3), (2, 4), (2, 5), (2, 6),
             (2, 7), (2, 8), (2, 9), (2, 10), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (3, 10), (4, 5), (4, 6),
             (4, 7), (4, 8), (4, 9), (4, 10), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10), (6, 7), (6, 8), (6, 9), (6, 10),
             (7, 8), (7, 9), (7, 10), (8, 9), (8, 10), (9, 10)], targetsplit._possible_helical_pairs)
        self.assertListEqual(
            ['1_2', '1_3', '1_4', '1_5', '1_6', '1_7', '1_8', '1_9', '1_10', '2_3', '2_4', '2_5', '2_6', '2_7', '2_8',
             '2_9', '2_10', '3_4', '3_5', '3_6', '3_7', '3_8', '3_9', '3_10', '4_5', '4_6', '4_7', '4_8', '4_9', '4_10',
             '5_6', '5_7', '5_8', '5_9', '5_10', '6_7', '6_8', '6_9', '6_10', '7_8', '7_9', '7_10', '8_9', '8_10',
             '9_10'], [x.id for x in targetsplit.subtargets])
        self.assertListEqual(['2_6', '1_7', '4_9', '3_5', '2_7', '3_4', '3_8', '4_10', '4_5', '8_10', '4_8', '9_10'],
                             [x.id for x in targetsplit.ranked_subtargets])
        self.assertTupleEqual(([42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62],
                               [164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
                                181, 182, 183, 184]), targetsplit.get_helical_pair((2, 6)))
        self.assertListEqual([100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
                              118, 119, 120], targetsplit._helix_range(4))
        self.assertListEqual(
            [(12, 36), (16, 40), (17, 40), (13, 39), (13, 40), (13, 36), (16, 37), (9, 36), (16, 41), (12, 33),
             (14, 40)], [x.id for x in targetsplit.get_interhelical_contacts((2, 6))])
