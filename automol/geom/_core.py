""" core library defining the geom data structure
"""


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
