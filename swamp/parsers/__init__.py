"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements classes and methods parse files of interest.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version
import os

__version__ = version.__version__

if 'THIS_IS_READTHEDOCS' not in os.environ and "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")


def AleigenParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.aleigenparser.AleigenParser` instance"""
    from swamp.parsers.aleigenparser import AleigenParser

    return AleigenParser(*args, **kwargs)


def MapAlignParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.mapalignparser.MapAlignParser` instance"""
    from swamp.parsers.mapalignparser import MapAlignParser

    return MapAlignParser(*args, **kwargs)


def Parser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.parser.Parser` instance"""
    from swamp.parsers.parser import Parser

    return Parser(*args, **kwargs)


def PdbtmXmlParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.pdbtmxmlparser.PdbtmXmlParser` instance"""
    from swamp.parsers.pdbtmxmlparser import PdbtmXmlParser

    return PdbtmXmlParser(*args, **kwargs)


def PhaserParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.phaser.PhaserParser` instance"""
    from swamp.parsers.phaserparser import PhaserParser

    return PhaserParser(*args, **kwargs)


def PhenixParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.phenixparser.PhenixParser` instance"""
    from swamp.parsers.phenixparser import PhenixParser

    return PhenixParser(*args, **kwargs)


def RefmacParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.refmacparser.RefmacParser` instance"""
    from swamp.parsers.refmacparser import RefmacParser

    return RefmacParser(*args, **kwargs)


def ShelxeParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.shelxeparser.ShelxeParser` instance"""
    from swamp.parsers.shelxeparser import ShelxeParser

    return ShelxeParser(*args, **kwargs)


def TopconsParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.topconsparser.TopconsParser` instance"""
    from swamp.parsers.topconsparser import TopconsParser

    return TopconsParser(*args, **kwargs)


def MtzParser(*args, **kwargs):
    """:py:obj:`~swamp.parsers.mtzparser.MtzParser` instance"""
    from swamp.parsers.mtzparser import MtzParser

    return MtzParser(*args, **kwargs)
