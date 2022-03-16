import openmesh as om

import sys
sys.path.append('../../')
from MeshUtility import ff_graph


def test_ff():
    mesh = om.read_trimesh('../data/sphere_cvt.obj')
    G = ff_graph(mesh)
    print(mesh.n_faces())
    print(G.number_of_nodes())
    # print([x for x in G])


if __name__ == '__main__':
    test_ff()
