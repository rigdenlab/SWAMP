#!/usr/bin/env ccp4-python

"""Script with the framework for SWAMP unittesting
Credits to Felix Simkovic for SIMBAD/scripts/run_tests.py
"""

from __future__ import print_function
import glob
import logging
import os
import sys
from unittest import TestLoader, TextTestRunner, TestSuite

SWAMP_DIR = os.path.join(os.path.dirname(__file__), "swamp")
PACKAGES = ["parsers", "utils"]


class SWAMPUnittestFramework(object):
    """Framework to run SWAMP unittesting"""

    def run(self, buffer=False, pattern="test*.py", verbosity=2):
        """Main routine for running the test cases"""

        tests = self.find_tests(SWAMP_DIR, pattern=pattern)
        if int(tests.countTestCases()) <= 0:
            msg = 'Could not find any tests to run in directory: {0}'.format(SWAMP_DIR) + os.linesep
            sys.stderr.write(msg)
            sys.exit(1)
        logging.disable(logging.CRITICAL)
        result = TextTestRunner(verbosity=verbosity, buffer=buffer).run(tests)
        logging.disable(logging.NOTSET)
        if result.wasSuccessful():
            exit(0)
        else:
            exit(1)

    def find_tests(self, directory, pattern="test*.py"):
        """Load a unittest test suite"""

        search_pattern = os.path.join(directory, "*")
        cases = [os.path.basename(folder) for folder in glob.iglob(search_pattern)
                 if os.path.isdir(folder) and os.path.basename(folder) in PACKAGES]

        return self._load_suite(cases, pattern, directory)

    def _load_suite(self, cases, pattern, directory):
        suite = TestSuite()
        for case in cases:
            path = os.path.join(directory, case, "tests")
            try:
                _suite = TestLoader().discover(path, pattern=pattern, top_level_dir=directory)
                suite.addTests(_suite)
                del _suite
            except ImportError:
                print("*** not a package: {0} ***".format(path))
        return suite


if __name__ == "__main__":

    # Mock CCP4 directories
    ccp4_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    ccp4_scr = os.path.abspath(os.path.join(os.path.dirname(__file__), "CCP4_SCR"))
    os.environ['CCP4'] = ccp4_root
    os.environ['CCP4_SCR'] = ccp4_scr
    if not os.path.isdir(ccp4_scr):
        os.mkdir(ccp4_scr)

    test = SWAMPUnittestFramework()
    test.run()
