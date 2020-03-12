import os
import unittest
import joblib
from swamp.utils import remove, create_tempfile, touch
from swamp.search.searchjob import SearchJob
from swamp.wrappers import MapAlign

CONPRED_DUMMY = """PFRMAT RR
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

LIBRARY = os.path.join(os.environ['CCP4_SCR'], 'library')


class SearchJobTestCase(unittest.TestCase):

    def test_1(self):
        fname = create_tempfile(CONPRED_DUMMY)
        os.mkdir(os.path.join(os.environ['CCP4_SCR'], 'library'))
        touch(os.path.join(LIBRARY, 'dummy_1.psicov'))
        touch(os.path.join(LIBRARY, 'dummy_2.psicov'))
        touch(os.path.join(LIBRARY, 'dummy_3.psicov'))
        self.addCleanup(remove, LIBRARY)
        self.addCleanup(remove, fname)
        search = SearchJob(id='test', workdir=os.path.join(os.environ['CCP4_SCR'], 'search_test'), query=fname,
                           template_library=LIBRARY, con_format="psicov", library_format="pdb", algorithm='mapalign',
                           python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'))
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'search_test'))
        self.assertTrue(os.path.isdir(os.path.join(os.environ['CCP4_SCR'], 'search_test')))
        self.assertEqual("""**********************************************************************
********************        SWAMP - SCANJOB         ******************
**********************************************************************

""", search.search_header)
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'search_test', 'tmp_query.eigenvct'),
                         search.tmp_eigen_query)
        self.assertIsNone(search._pdbfile_template('test'))
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'search_test',
                                      os.path.basename(os.path.basename('4ryr_5B_3A', ).split(".")[0])),
                         search._job_dir_template('4ryr_5B_3A', ))
        self.assertEqual(os.path.join(os.environ['CCP4_SCR'], 'search_test', 'search_test_results.pckl'),
                         search.pickle_fname)
        self.assertListEqual([os.path.join(LIBRARY, 'dummy_1.psicov'), os.path.join(LIBRARY, 'dummy_2.psicov'),
                              os.path.join(LIBRARY, 'dummy_3.psicov')], sorted(search.template_list))

    def test_2(self):
        fname = create_tempfile(CONPRED_DUMMY)
        os.mkdir(os.path.join(os.environ['CCP4_SCR'], 'library'))
        touch(os.path.join(LIBRARY, '4ryr_3A_5B.psicov'))
        touch(os.path.join(LIBRARY, '3txt_6A_11B.psicov'))
        touch(os.path.join(LIBRARY, '1hvv_4A_57.psicov'))
        self.addCleanup(remove, LIBRARY)
        self.addCleanup(remove, fname)
        search = SearchJob(id='test', workdir=os.path.join(os.environ['CCP4_SCR'], 'search_test'), query=fname,
                           template_subset=('4ryr_5B_3A',), template_library=LIBRARY,
                           con_format="mapalign", library_format="mapalign", algorithm='mapalign',
                           python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'))
        self.addCleanup(remove, os.path.join(os.environ['CCP4_SCR'], 'search_test'))
        self.assertListEqual([os.path.join(LIBRARY, '4ryr_3A_5B.psicov')], search.template_list)
        self.assertIsInstance(MapAlign, type(search._alignment_wrapper))
        self.assertEqual("""cd %s
%s << EOF
from swamp.search import SearchJob
job=SearchJob(id="test", workdir="%s", query="%s", template_library="%s", con_format="mapalign", library_format="mapalign", template_subset=("4ryr_5B_3A"), algorithm="mapalign")
job.run()
job.store_pickle()
EOF
""" % (os.path.join(os.environ['CCP4_SCR'], 'search_test'), os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'),
       os.path.join(os.environ['CCP4_SCR'], 'search_test'), fname, LIBRARY), search._python_script)

        self.assertEqual(1, len(search.joblist))
        alignment_info = search._get_alignment_info(search.template_list[0])
        del alignment_info['logger']
        self.assertDictEqual(
            {'map_b': '%s/4ryr_3A_5B.psicov' % LIBRARY, 'map_a': fname, 'format_a': 'mapalign',
             'format_b': 'mapalign', 'workdir': os.path.join(os.environ['CCP4_SCR'], 'search_test')}, alignment_info)

    def test_3(self):
        fname = create_tempfile(CONPRED_DUMMY)
        os.mkdir(os.path.join(os.environ['CCP4_SCR'], 'library'))
        touch(os.path.join(LIBRARY, 'dummy_1.psicov'))
        touch(os.path.join(LIBRARY, 'dummy_2.psicov'))
        touch(os.path.join(LIBRARY, 'dummy_3.psicov'))
        self.addCleanup(remove, LIBRARY)
        self.addCleanup(remove, fname)
        search = SearchJob(id='test', workdir=os.path.join(os.environ['CCP4_SCR'], 'search_test'), query=fname,
                           template_library=LIBRARY, con_format="psicov", library_format="pdb", algorithm='mapalign',
                           python_interpreter=os.path.join(os.environ['CCP4'], 'bin', 'ccp4-python'))

        self.assertListEqual([], search.results)
        search.results = ['test_1', 'test_2']
        search.store_pickle()
        self.assertTrue(os.path.isfile(os.path.join(os.environ['CCP4_SCR'], 'search_test', 'search_test_results.pckl')))
        self.addCleanup(remove, search.pickle_fname)
        test_results = joblib.load(os.path.join(os.environ['CCP4_SCR'], 'search_test', 'search_test_results.pckl'))
        self.assertListEqual(['test_1', 'test_2'], test_results)
