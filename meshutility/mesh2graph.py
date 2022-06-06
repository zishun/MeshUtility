"""
=============
MESH TO GRAPH
=============

Construct mesh vert-vert graph and facet-facet graph.

by Zishun Liu <liuzishun@gmail.com>
Feb 15, 2021
"""

import openmesh as om
import networkx as nx

__all__ = ['ff_graph', 'vv_graph']


def ff_graph(mesh):
    G_face = nx.empty_graph(mesh.n_faces())
    indices = mesh.ef_indices()
    G_face.add_edges_from(indices)
    if G_face.has_node(-1):
        G_face.remove_node(-1)
    return G_face


def vv_graph(mesh, edge_length_weight=False):
    G = nx.empty_graph(mesh.n_vertices())
    indices = mesh.ev_indices()
    if not edge_length_weight:
        G.add_edges_from(indices)
    else:
        pts = mesh.points()
        lengths = np.linalg.norm(pts[indices[:,0]]-pts[indices[:,1]], axis=1)
        edges = [(e[0], e[1], l) for e, l in zip(indices, lengths)]
        G.add_weighted_edges_from(edges)
    return G


if __name__ == '__main__':
    pass
