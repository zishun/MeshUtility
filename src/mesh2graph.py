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


def vv_graph(mesh):
    G = nx.empty_graph(mesh.n_vertices())
    indices = mesh.ev_indices()
    G.add_edges_from(indices)
    return G


if __name__ == '__main__':
    pass
