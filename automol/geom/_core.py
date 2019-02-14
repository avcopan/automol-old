""" core library defining the geom data structure
"""
from .. import atom as _atom


def from_symbols_and_coordinates(syms, xyzs):
    """ geometry data structure from symbols and coordinates
    """
    syms = list(map(str.upper, syms))
    assert all(sym in _atom.SYMBOLS for sym in syms)
    geo = tuple(zip(syms, xyzs))
    return geo


def symbols(geo):
    """ a dictionary of atomic symbols by index
    """
    syms, _ = zip(*geo)
    return syms


def coordinates(geo):
    """ a dictionary of coordinates by index
    """
    _, xyzs = zip(*geo)
    return xyzs
