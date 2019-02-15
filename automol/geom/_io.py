""" I/O for cartesian geometries
"""
import numpy
import autoparse.pattern as app
import autoparse.find as apf
from .._cnst.geom import from_data
from .. import _units
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates

ATOM_SYMBOL_PATTERN = app.LETTER + app.maybe(app.LETTER)


def string(geo, to_angstroms=True):
    """ write the cartesian geometry as a string
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo)
    xyzs = xyzs if not to_angstroms else numpy.multiply(xyzs, _units.BOHR2ANG)

    geo_str = '\n'.join('{:2s} {:20.12f} {:20.12f} {:20.12f}'.format(sym, *xyz)
                        for sym, xyz in zip(syms, xyzs))
    return geo_str


def from_string(geo_str, angstroms=True, strict=True):
    """ read a cartesian geometry from a string
    """
    sym_capturing_pattern = app.LINESPACES.join(
        [app.capturing(ATOM_SYMBOL_PATTERN)] + [app.FLOAT] * 3)
    xyz_capturing_pattern = app.LINESPACES.join(
        [ATOM_SYMBOL_PATTERN] + [app.capturing(app.FLOAT)] * 3)

    if strict:
        line_pattern = app.maybe(app.LINESPACES).join(
            [app.LINE_START, ATOM_SYMBOL_PATTERN] + [app.FLOAT] * 3 +
            [app.LINE_END])
        lines = geo_str.splitlines()
        assert all(apf.has_match(line_pattern, line) for line in lines)

    syms = apf.all_captures(sym_capturing_pattern, geo_str)
    xyz_strs_lst = apf.all_captures(xyz_capturing_pattern, geo_str)
    xyzs = tuple(tuple(map(float, xyz_strs)) for xyz_strs in xyz_strs_lst)
    geo = from_data(syms, xyzs, angstroms=angstroms)
    return geo


def dxyz_string(geo, comment_line=''):
    """ write the cartesian geometry to a .xyz string
    """
    natms = len(_symbols(geo))
    assert not apf.has_match(app.NEWLINE, comment_line)
    geo_str = string(geo)
    dxyz_str = '{:d}\n{:s}\n{:s}'.format(natms, comment_line, geo_str)
    return dxyz_str


def from_dxyz_string(dxyz_str, with_comment_line=False):
    """ read a cartesian geometry from a .xyz string
    """
    lines = dxyz_str.splitlines()
    assert apf.has_match(app.UNSIGNED_INTEGER, lines[0])
    natms = int(lines[0])
    comment_line = lines[1]
    geo_str = '\n'.join(lines[2:natms+2])
    geo = from_string(geo_str, angstroms=True, strict=True)
    return geo if not with_comment_line else (geo, comment_line)
