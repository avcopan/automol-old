""" some math functions
"""
import numpy
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def almost_equal(geom1, geom2):
    """ are these geometries almost equal
    """
    ret = False
    if _symbols(geom1) == _symbols(geom2):
        ret = numpy.allclose(_coordinates(geom1), _coordinates(geom2))
    return ret
