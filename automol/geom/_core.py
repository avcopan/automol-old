""" core library defining the geom data structure
"""
from numbers import Real as _Float
from .. import atom as _atom


def from_symbols_and_coordinates(syms, xyzs):
    """ geometry data structure from symbols and coordinates
    """
    assert len(syms) == len(xyzs)
    geo = tuple(_atom_entry(sym, xyz) for sym, xyz in zip(syms, xyzs))
    return geo


def _atom_entry(sym, xyz):
    """ enforce the correct format for an atom entry
    """
    sym = _atom.standard_case(sym)
    assert sym in _atom.SYMBOLS
    assert len(xyz) == 3
    assert all(isinstance(comp, _Float) for comp in xyz)
    xyz = tuple(map(float, xyz))
    return (sym, xyz)


def symbols(geo):
    """ a dictionary of atomic symbols by index
    """
    if geo:
        syms, _ = zip(*geo)
    else:
        syms = ()
    return syms


def coordinates(geo):
    """ a dictionary of coordinates by index
    """
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    return xyzs
