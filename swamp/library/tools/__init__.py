"""This is SWAMP: Solving structures With Alpha Membrane Pairs

This module implements classes and methods to create and store the data contained in the library of helical pairs.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version

__version__ = version.__version__

import os
import gzip

if "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")


def compress(fname, out=None):
    """Compress a text file into .gz

    :param fname: the file name to be compressed
    :type fname: str
    :param out: specify an output file name, otherwise default is fname.gz
    :type out: str, None
    :return compressed file name
    :rtype str
    """

    if out is None:
        out = '%s.gz' % fname

    with open(fname, 'rb') as f_in, gzip.open(out, 'wb') as f_out:
        data = f_in.read()
        bindata = bytearray(data)
        f_out.write(bindata)

    return out


def decompress(fname, out=None):
    """Decompress a .gz file into text file

    :param fname: the file name to be decompressed
    :type fname: str
    :param out: specify an output file name, otherwise default is fname without .gz
    :type out: str, None
    :return decompressed file name
    :rtype str
    """

    if out is None:
        out = fname.replace('.gz', '')

    with open(out, "wb") as f_out, gzip.open(fname, "rb") as f_in:
        bindata = f_in.read()
        f_out.write(bindata)

    return out
