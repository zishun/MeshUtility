"""
==========
SPLIT MESH
==========

Split mesh.

by Zishun Liu <liuzishun@gmail.com>
Feb 13, 2021
"""

import numpy as np
import openmesh as om

__all__ = ['split_mesh']


def split_mesh(V, F, edges, ratios):
    """
    Split mesh.

    Parameters
    ----------
    V : numpy.array
    F : numpy.array
    edges : numpy.array
    ratios : numpy.array
        p0 * (1-r)+ p1 * r
    Returns
    -------
    mesh : openmesh.TriMesh
    index : list (int)
    """

    mesh = om.TriMesh(V, F)

    points = mesh.points()
    pts0 = points[edges[:, 0], :]
    pts1 = points[edges[:, 1], :]
    curve_pts = np.multiply(pts0, 1.-ratios[:, np.newaxis]) + \
        np.multiply(pts1, ratios[:, np.newaxis])

    new_index = []
    for i in range(edges.shape[0]):
        if edges[i, 0] == edges[i, 1]:
            new_index.append(edges[i, 0])
        else:
            p = curve_pts[i, :]
            heh = mesh.find_halfedge(mesh.vertex_handle(edges[i, 0]),
                                     mesh.vertex_handle(edges[i, 1]))
            eh = mesh.edge_handle(heh)
            vh = mesh.add_vertex(p)
            mesh.split_edge(eh, vh)
            new_index.append(vh.idx())

    return mesh, new_index


if __name__ == '__main__':
    pass
