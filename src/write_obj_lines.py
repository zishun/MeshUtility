"""
===============
WRITE OBJ LINES
===============

Author: Zishun Liu <liuzishun@gmail.com>
Date: Feb 16, 2021
"""


def write_obj_lines(fn, verts, lines):
    """
    Write polylines to an obj file.
    Note: Each line in the file only contains one segment so that Meshlab
    opens it properly.

    Parameters
    ----------
    fn : str
        path to the output obj file
    verts : numpy.array
        vertex coordinates (n x 3)
    lines : list of numpy.array
        list of curves [arr_0, arr_1, ...]
        The i-th curve, arr_i, is a vector of vertex indices of vertices on it,
        0-based.
    Returns
    -------
    None
    """
    with open(fn, 'w') as f:
        for v in verts:
            f.write('v %.8g %.8g %.8g\n' % (v[0], v[1], v[2]))
        for line in lines:
            for i in range(len(line)-1):
                f.write('l %d %d\n' % (line[i]+1, line[i+1]+1))
