"""This is SWAMP: Solving Structures With Alpha Membrane Pairs

This module implements wrapper classes around programs used within SWAMP mr module.
"""

__author__ = "Filomeno Sanchez Rodriguez"
__credits__ = "Daniel Rigden, & Ronan Keegan"
__email__ = "filomeno.sanchez-rodriguez@liv.ac.uk"

from swamp import version
import os

__version__ = version.__version__

if 'THIS_IS_READTHEDOCS' not in os.environ and "CCP4" not in os.environ:
    raise RuntimeError("Cannot find CCP4 root directory")


def AlEigen(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.aleigen.AlEigen` instance"""
    from swamp.wrappers.aleigen import AlEigen

    return AlEigen(*args, **kwargs)


def MapAlign(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.mapalign.MapAlign` instance"""
    from swamp.wrappers.mapalign import MapAlign

    return MapAlign(*args, **kwargs)


def Wrapper(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.wrapper.Wrapper` instance"""
    from swamp.wrappers.wrapper import Wrapper

    return Wrapper(*args, **kwargs)


def Phaser(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.wphaser.Phaser` instance"""
    from swamp.wrappers.wphaser import Phaser

    return Phaser(*args, **kwargs)


def PhenixCC(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.PhenixCC` instance"""
    from swamp.wrappers.phenixcc import PhenixCC

    return PhenixCC(*args, **kwargs)


def wRefmac(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.wrefmac.wRefmac` instance"""
    from swamp.wrappers.wrefmac import wRefmac

    return wRefmac(*args, **kwargs)


def Shelxe(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.shelxe.Shelxe` instance"""
    from swamp.wrappers.shelxe import Shelxe

    return Shelxe(*args, **kwargs)


def Mtz2Various(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.mtz2various.Mtz2Various` instance"""
    from swamp.wrappers.mtz2various import Mtz2Various

    return Mtz2Various(*args, **kwargs)

def Gesamt(*args, **kwargs):
    """:py:obj:`~swamp.wrappers.gesamt.Gesamt` instance"""
    from swamp.wrappers.gesamt import Gesamt

    return Gesamt(*args, **kwargs)
