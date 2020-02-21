import os
import unittest
from swamp.wrappers.mtz2various import Mtz2Various


class MyTestCase(unittest.TestCase):

    def test_1(self):
        mtz2various = Mtz2Various(mtzin='/empty/path/fname.mtz', hklout='/empty/path/fname.hkl', workdir='/empty/path')
        self.assertListEqual(mtz2various.cmd,
                             ["mtz2various", 'HKLIN', '/empty/path/fname.mtz', 'HKLOUT', '/empty/path/fname.hkl'])


if __name__ == '__main__':
    unittest.main()
