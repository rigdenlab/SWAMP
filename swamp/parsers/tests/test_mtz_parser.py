import unittest
from swamp.parsers.mtzparser import MtzParser


class MockMtzParser(MtzParser):
    """A class to mock :py:obj:`~swmap.parsers.mtzparser.MtzParser` for testing purposes"""

    def read_reflections(self):
        """Override :py:func:`~swamp.parsers.mtzparser.Mtzparser.read_reflections` for testing purposes"""
        self._all_labels = []

    @property
    def all_labels(self):
        """Setter for :py:attr:`~swamp.parsers.mtzparser.MtzParser.all_labels`"""
        return self._all_labels

    @all_labels.setter
    def all_labels(self, value):
        self._all_labels = value


class MtzParserTestCase(unittest.TestCase):

    def test_1(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        parser.all_labels = ['SIGI(-)', 'FP', 'F(+)', 'fp(-)', 'sigF(+)', 'sigi(+)', 'SIGFP(-)',
                             'dummy-label', 'Free', 'SIGFP', 'sigi', 'i', 'I(+)', 'i(-)']
        parser.error = False
        parser.parse()
        self.assertTupleEqual(parser.summary, (b'FP', b'SIGFP', b'i', b'sigi', b'Free', b'F(+)', b'sigF(+)', b'I(+)',
                                               b'sigi(+)', b'fp(-)', b'SIGFP(-)', b'i(-)', b'SIGI(-)'))

    def test_2(self):
        parser = MockMtzParser('/empty/path/fname.mtz')
        parser.all_labels = ['dummy-label']
        self.assertTrue(parser.error)
        parser.error = False
        parser.parse()
        self.assertTrue(parser.error)
