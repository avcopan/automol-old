""" geom constructor
"""
from numbers import Real as _Real
import numpy
from .. import atom as _atom
from .. import _units


def from_data(syms, xyzs, angstroms=False):
    """ geometry data structure from symbols and coordinates
    """
    assert len(syms) == len(xyzs)
    xyzs = xyzs if not angstroms else numpy.multiply(xyzs, _units.ANG2BOHR)
    geo = tuple(_atom_entry(sym, xyz) for sym, xyz in zip(syms, xyzs))
    return geo


def _atom_entry(sym, xyz):
    """ enforce the correct format for an atom entry
    """
    sym = _atom.standard_case(sym)
    assert sym in _atom.SYMBOLS
    assert len(xyz) == 3
    assert all(isinstance(comp, _Real) for comp in xyz)
    xyz = tuple(map(float, xyz))
    return (sym, xyz)
