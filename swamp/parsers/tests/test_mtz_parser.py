import unittest
from swamp.parsers.mtzparser import MtzParser


class MockGemmiMtz(object):
    """A class to mock :py:obj:gemmi.Mtz for testing purposes

    :param list columns: a list containing instances of :py:obj:`~swmap.parsers.tests.test_mtz_parser.MockGemmiMtzColumn`
    """

    def __init__(self, columns):
        self.columns = columns


class MockGemmiMtzColumn(object):
    """A class to mock :py:obj:gemmi.Mtz.Column for testing purposes

    :param str label: label of the mocked MTZ column
    :param str type: type of the mocked MTZ column
    """

    def __init__(self, label, type):
        self.label = label
        self.type = type


class MockMtzParser(MtzParser):
    """A class to mock :py:obj:`~swmap.parsers.mtzparser.MtzParser` for testing purposes"""

    def read_reflections(self):
        """Override :py:func:`~swamp.parsers.mtzparser.Mtzparser.read_reflections` for testing purposes"""
        pass


class MtzParserTestCase(unittest.TestCase):

    def test_1(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        mock_columns = [MockGemmiMtzColumn("SIGI(-)", "M"),
                        MockGemmiMtzColumn("FP", "F"),
                        MockGemmiMtzColumn("F(+)", "G"),
                        MockGemmiMtzColumn("fp(-)", "G"),
                        MockGemmiMtzColumn("sigF(+)", "L"),
                        MockGemmiMtzColumn("sigi(+)", "M"),
                        MockGemmiMtzColumn("SIGFP(-)", "L"),
                        MockGemmiMtzColumn("dummy-label", "E"),
                        MockGemmiMtzColumn("Free", "I"),
                        MockGemmiMtzColumn("SIGFP", "Q"),
                        MockGemmiMtzColumn("sigi", "Q"),
                        MockGemmiMtzColumn("i", "J"),
                        MockGemmiMtzColumn("I(+)", "K"),
                        MockGemmiMtzColumn("i(-)", "K")]
        parser.reflection_file = MockGemmiMtz(mock_columns)
        parser.error = False
        parser.parse()
        self.assertTupleEqual(parser.summary, (b'FP', b'SIGFP', b'i', b'sigi', b'Free', b'F(+)', b'sigF(+)', b'I(+)',
                                               b'sigi(+)', b'fp(-)', b'SIGFP(-)', b'i(-)', b'SIGI(-)'))
        self.assertEqual(parser.fname, '/empty/path/fname.mtz')

    def test_2(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        mock_columns = [MockGemmiMtzColumn("dummy-label", "E")]
        parser.reflection_file = MockGemmiMtz(mock_columns)
        self.assertTrue(parser.error)
        parser.error = False
        parser.parse()
        self.assertTrue(parser.error)

    def test_3(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        mock_columns = [MockGemmiMtzColumn("H", "H"),
                        MockGemmiMtzColumn("K", "H"),
                        MockGemmiMtzColumn("L", "H"),
                        MockGemmiMtzColumn("FREE", "I"),
                        MockGemmiMtzColumn("FP", "F"),
                        MockGemmiMtzColumn("SIGFP", "Q"),
                        MockGemmiMtzColumn("FC", "F"),
                        MockGemmiMtzColumn("PHIC", "P"),
                        MockGemmiMtzColumn("FC_ALL", "F"),
                        MockGemmiMtzColumn("PHIC_ALL", "P"),
                        MockGemmiMtzColumn("FWT", "F"),
                        MockGemmiMtzColumn("PHWT", "P"),
                        MockGemmiMtzColumn("DELFWT", "F"),
                        MockGemmiMtzColumn("PHDELWT", "P"),
                        MockGemmiMtzColumn("FOM", "W"),
                        MockGemmiMtzColumn("FC_ALL_LS", "F"),
                        MockGemmiMtzColumn("PHIC_ALL_LS", "P")]
        parser.reflection_file = MockGemmiMtz(mock_columns)
        parser.error = False
        parser.parse()
        print(parser.summary)
        self.assertTupleEqual(parser.summary, (b'FP', b'SIGFP', None, None, b'FREE', None, None, None,
                                               None, None, None, None, None))
        self.assertEqual(parser.fname, '/empty/path/fname.mtz')

    def test_4(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        mock_columns = [MockGemmiMtzColumn("H", "H"),
                        MockGemmiMtzColumn("K", "H"),
                        MockGemmiMtzColumn("L", "H"),
                        MockGemmiMtzColumn("FreeR_flag", "I"),
                        MockGemmiMtzColumn("F_XDSdataset", "F"),
                        MockGemmiMtzColumn("SIGF_XDSdataset", "Q"),
                        MockGemmiMtzColumn("DANO_XDSdataset", "D"),
                        MockGemmiMtzColumn("SIGDANO_XDSdataset", "Q"),
                        MockGemmiMtzColumn("F_XDSdataset(+)", "G"),
                        MockGemmiMtzColumn("SIGF_XDSdataset(+)", "L"),
                        MockGemmiMtzColumn("F_XDSdataset(-)", "G"),
                        MockGemmiMtzColumn("SIGF_XDSdataset(-)", "L"),
                        MockGemmiMtzColumn("ISYM_XDSdataset", "Y"),
                        MockGemmiMtzColumn("IMEAN_XDSdataset", "J"),
                        MockGemmiMtzColumn("SIGIMEAN_XDSdataset", "Q"),
                        MockGemmiMtzColumn("I_XDSdataset(+)", "K"),
                        MockGemmiMtzColumn("SIGI_XDSdataset(+)", "M"),
                        MockGemmiMtzColumn("I_XDSdataset(-)", "K"),
                        MockGemmiMtzColumn("SIGI_XDSdataset(-)", "M")]
        parser.reflection_file = MockGemmiMtz(mock_columns)
        parser.error = False
        parser.parse()

        self.assertTupleEqual(parser.summary, (b'F_XDSdataset', b'SIGF_XDSdataset', b'IMEAN_XDSdataset',
                                               b'SIGIMEAN_XDSdataset', b'FreeR_flag', b'F_XDSdataset(+)',
                                               b'SIGF_XDSdataset(+)', b'I_XDSdataset(+)', b'SIGI_XDSdataset(+)',
                                               b'F_XDSdataset(-)', b'SIGF_XDSdataset(-)', b'I_XDSdataset(-)',
                                               b'SIGI_XDSdataset(-)'))
        self.assertEqual(parser.fname, '/empty/path/fname.mtz')
