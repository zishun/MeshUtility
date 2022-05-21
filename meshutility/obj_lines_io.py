"""
===============
WRITE OBJ LINES
===============

Author: Zishun Liu <liuzishun@gmail.com>
Date: Feb 16, 2021
"""
import numpy as np


__all__ = ['write_obj_lines', 'read_obj_lines']


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


def read_obj_lines(fn):
    """
    Read polylines from an obj file.
    
    Parameters
    ----------
    fn : str
        path to the output obj file
    Returns
    -------
    verts : numpy.array
        vertex coordinates (n x 3)
    lines : numpy.array (m x 2)
        edge segments
    """
    verts = []
    segs = []
    with open(fn, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        s = lines[i].split()
        if (len(s) == 0):
            continue
        if (s[0] == 'v'):
            verts.append([float(s[1]), float(s[2]), float(s[3])])
        elif (s[0] == 'l'):
            L = list(map(lambda x: int(x),  s[1:]))
            segs.extend([[L[i], L[i+1]] for i in range(len(L)-1)])

    verts_mat = np.array(verts)
    segs_mat = np.array(segs)-1  # 1-indexed -> 0-indexed
    return verts_mat, segs_mat
