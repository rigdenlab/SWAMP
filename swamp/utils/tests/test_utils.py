import gemmi
import unittest
from swamp.utils import *
from conkit.core import ContactMap, Contact


class UtilsTestCase(unittest.TestCase):

    def test_contactmap_operations(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence("TEST", "ACDEFGHK")
        inverted = invert_contactmap(contact_map)
        self.assertListEqual([x.id for x in inverted], [(8, 5), (7, 5), (4, 1), (6, 3)])

    def test_gemmi_hierarchy_operations(self):
        pdb_content = """CRYST1   73.330   73.330  163.520  90.00  90.00  90.00 P 41 2 2      8          
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
"""
        to_extract = [1, 2, 3, 4, 5, 9]
        hierarchy = gemmi.read_pdb_string(pdb_content)
        renumber_hierarchy(hierarchy)
        self.assertListEqual(list(range(1, 12)), [x.seqid.num for x in hierarchy[0][0]])
        new_hierarchy = extract_hierarchy(hierarchy, to_extract)
        self.assertListEqual([x.seqid.num for x in new_hierarchy[0][0]], to_extract)
        inverted_hierarchy = invert_hiearchy(new_hierarchy)
        self.assertListEqual([x.name for x in inverted_hierarchy[0][0]],
                             list(reversed([y.name for y in new_hierarchy[0][0]])))
        merged_hierarchy = merge_into_ensemble((inverted_hierarchy, new_hierarchy))
        self.assertEqual(len(merged_hierarchy), 2)
        self.assertListEqual([x.name for x in merged_hierarchy[0][0]],
                             list([y.name for y in inverted_hierarchy[0][0]]))
        self.assertListEqual([x.name for x in merged_hierarchy[1][0]],
                             list([y.name for y in new_hierarchy[0][0]]))
        models = split_ensemble_into_models(merged_hierarchy)
        self.assertListEqual([x.name for x in models[0][0][0]],
                             list([y.name for y in inverted_hierarchy[0][0]]))
        self.assertListEqual([x.name for x in models[1][0][0]],
                             list([y.name for y in new_hierarchy[0][0]]))
        merged = merge_hierarchies((inverted_hierarchy, new_hierarchy))
        self.assertListEqual([x.seqid.num for x in merged[0][0]],
                             [x.seqid.num for x in inverted_hierarchy[0][0]]
                             + [x.seqid.num for x in new_hierarchy[0][0]])

        helices = ((1, 2, 3, 4), (5, 6, 7, 8))
        fragment_cmap = extract_fragment_cmap(hierarchy, helices)
        self.assertListEqual([x.id for x in fragment_cmap],
                             [(1, 5), (2, 5), (2, 6), (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (4, 8)])
