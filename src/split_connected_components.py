"""
==========================
SPLIT CONNECTED COMPONENTS
==========================

Split connected components.

by Zishun Liu <liuzishun@gmail.com>
May 10, 2021
"""

import numpy as np
import openmesh as om
import networkx as nx

from .mesh2graph import ff_graph

__all__ = ['split_connected_components']


def split_connected_components(V, F):
    """
    Split connected components.

    Parameters
    ----------
    V : numpy.array
    F : numpy.array
    Returns
    -------
    parts : list (openmesh.PolyMesh)
    """

    mesh = om.PolyMesh(V, F)
    G = ff_graph(mesh)

    i = 0
    parts = []
    for c in nx.connected_components(G):
        mask = np.zeros((mesh.n_faces(),), dtype=bool)
        for x in c:
            mask[x] = True
        part = om.PolyMesh(V, F[mask])
        for vh in part.vertices():
            if not any(True for _ in part.vf(vh)):
            # if sum([1 for _ in part.vf(vh)]) == 0:
                part.delete_vertex(vh)
        part.garbage_collection()
        parts.append(part)
        i += 1

    return parts


if __name__ == '__main__':
    pass
