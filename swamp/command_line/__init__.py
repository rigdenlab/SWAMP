"""This is SWAMP: Solving structure With Alpha Membrane Pairs

This module contains python scripts that can be used top perform the most commonly used tasks from the command line
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version
import os

__version__ = version.__version__

if 'THIS_IS_READTHEDOCS' not in os.environ and "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")


def check_file_exists(input_path):
    """Check if a given path corresponds with an existing file

    :param input_path: location of the file to be tested
    :type input_path: str, None
    :returns the absolute path of the file if it exists, None if the input is None
    :rtype str, None
    :raises IOError if the file doesn't exist
    """

    if input_path is None:
        return None
    if os.path.exists(os.path.abspath(input_path)):
        return os.path.abspath(input_path)
    else:
        raise IOError("%s file/directory does not exist!" % input_path)
