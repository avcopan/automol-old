""" graph conversion helpers
"""
import itertools
import numpy
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def connectivity_graph(geo):
    """ the connectivity graph
    """
    # using the same cut-offs as x2z:
    xy_bond_max = 3.5 / 1.8897259886
    xh_bond_max = 2.5 / 1.8897259886

    atm_sym_dct = dict(enumerate(_symbols(geo)))
    atm_xyz_dct = dict(enumerate(_coordinates(geo)))

    def _are_bonded(atm_key1, atm_key2):
        atm_sym1 = atm_sym_dct[atm_key1]
        atm_sym2 = atm_sym_dct[atm_key2]
        atm_xyz1 = atm_xyz_dct[atm_key1]
        atm_xyz2 = atm_xyz_dct[atm_key2]
        dist = numpy.linalg.norm(numpy.subtract(atm_xyz1, atm_xyz2))

        ret = False
        if 'H' in (atm_sym1, atm_sym2):
            ret = (dist < xh_bond_max)
        else:
            ret = (dist < xy_bond_max)

        return ret

    atm_dct = {atm_key: (atm_sym, 0, None)
               for atm_key, atm_sym in atm_sym_dct.items()}
    bnd_dct = {frozenset({atm_key1, atm_key2}): (1, None)
               for (atm_key1, atm_key2)
               in itertools.combinations(atm_dct.keys(), r=2)
               if _are_bonded(atm_key1, atm_key2)}
    cgr = (atm_dct, bnd_dct)
    return cgr
