""" test the automol.geom module
"""
from automol import geom

C2H2CLF_ICH = 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'
C2H2CLF_STE_ICH = 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1+'
C2H2CLF_GEO = (('F', (1.584822920001, -0.748486564300, -0.4271224303432)),
               ('C', (0.619219854789, 0.190165523008, -0.2716389279610)),
               ('C', (-0.635730620967, -0.1839138961594, -0.1803643636082)),
               ('Cl', (-1.602333181611, 0.736677675476, -0.02605091648865)),
               ('H', (0.916321356258, 1.229945559249, -0.2271265738271)),
               ('H', (-0.882300328471, -1.224388297273, -0.229635969682)))
C2H2CLF_CGR = ({0: ('F', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
                3: ('Cl', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
                frozenset({1, 4}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({2, 5}): (1, None)})
C2H2CLF_DXYZ_STR = """6
charge: 0, mult: 1
F        1.584822920001      -0.748486564300      -0.427122430343
C        0.619219854789       0.190165523008      -0.271638927961
C       -0.635730620967      -0.183913896159      -0.180364363608
Cl      -1.602333181611       0.736677675476      -0.026050916489
H        0.916321356258       1.229945559249      -0.227126573827
H       -0.882300328471      -1.224388297273      -0.229635969682
"""


def test__from_symbols_and_coordinates():
    """ test geom.from_symbols_and_coordinates
    """
    assert C2H2CLF_GEO == geom.from_symbols_and_coordinates(
        syms=geom.symbols(C2H2CLF_GEO),
        xyzs=geom.coordinates(C2H2CLF_GEO),
    )


def test__connectivity_graph():
    """ test geom.connectivity_graph
    """
    assert geom.connectivity_graph(C2H2CLF_GEO) == C2H2CLF_CGR


def test__inchi():
    """ test geom.inchi
    """
    assert geom.inchi(C2H2CLF_GEO) == C2H2CLF_ICH


def test__stereo_inchi():
    """ test geom.inchi
    """
    assert geom.stereo_inchi(C2H2CLF_GEO) == C2H2CLF_STE_ICH


def test__from_string():
    """ test geom.from_string
    """
    assert geom.almost_equal(
        geom.from_string(geom.string(C2H2CLF_GEO)), C2H2CLF_GEO)


def test__from_dxyz_string():
    """ test geom.from_dxyz_string
    """
    assert geom.almost_equal(
        geom.from_dxyz_string(geom.dxyz_string(C2H2CLF_GEO)), C2H2CLF_GEO)

    geo, comment_line = geom.from_dxyz_string(C2H2CLF_DXYZ_STR,
                                              with_comment_line=True)
    dxyz_str = geom.dxyz_string(geo, comment_line=comment_line)
    assert dxyz_str.rstrip() == C2H2CLF_DXYZ_STR.rstrip()


if __name__ == '__main__':
    # test__connectivity_graph()
    # test__inchi()
    # test__stereo_inchi()
    # test__from_string()
    # test__from_dxyz_string()
    test__from_symbols_and_coordinates()
