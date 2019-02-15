""" core library defining the geom data structure
"""


def symbols(geo):
    """ atomic symbols
    """
    if geo:
        syms, _ = zip(*geo)
    else:
        syms = ()
    return syms


def coordinates(geo):
    """ atomic coordinates
    """
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    return xyzs
