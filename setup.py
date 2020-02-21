#!/usr/bin/env ccp4-python

"""Script to setup necessary dependencies for SWAMP and manage unit testing"""

from __future__ import print_function
import argparse
import glob
import logging
import os
import sys
import subprocess
from unittest import TestLoader, TextTestRunner, TestSuite

SWAMP_DIR = os.path.join(os.path.dirname(__file__), "swamp")
REQUIREMENTS = os.path.join(os.path.dirname(__file__), "requirements.txt")
PACKAGES = ["parsers", "utils", "wrappers"]


def parse_arguments():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='SWAMP-SETUP: Setup necessary dependencies'
                                                 ' for SWAMP and manage unit testing')
    parser.add_argument("--tests", action='store_true',
                        help='If set, run SWAMP unit testing')

    args = parser.parse_args()

    return args


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

    args = parse_arguments()

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pybind11'], stdout=open(os.devnull, 'wb'))
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-r', REQUIREMENTS], stdout=open(os.devnull, 'wb'))

    if args.tests:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'codecov'], stdout=open(os.devnull, 'wb'))
        # Mock CCP4 directories for Travis CI
        ccp4_root = os.path.abspath(os.path.join(os.path.dirname(__file__)))
        ccp4_scr = os.path.abspath(os.path.join(os.path.dirname(__file__), "CCP4_SCR"))
        os.environ['CCP4'] = ccp4_root
        os.environ['CCP4_SCR'] = ccp4_scr
        if not os.path.isdir(ccp4_scr):
            os.mkdir(ccp4_scr)

        test = SWAMPUnittestFramework()
        test.run()
