import openmesh as om

import sys
sys.path.append('../../')
from MeshUtility import pygeodesic, colormap_vertex_color, \
split_mesh


def test_geodesic():
    # run python sphere_cvt.py to generate the mesh
    mesh = om.read_trimesh('../data/sphere_cvt.obj')
    # find source
    u, v = 0, 20
    path_edge, path_ratio = pygeodesic.find_path(mesh.points(),
                                                 mesh.face_vertex_indices(),
                                                 u, v)
    #source = split_along_curve(mesh, path_edge, path_ratio)
    mesh, source = split_mesh(mesh.points(), mesh.face_vertex_indices(),
                              path_edge, path_ratio)

    # compute geodesic distance field
    field = pygeodesic.distance_field(mesh.points(),
                                      mesh.face_vertex_indices(),
                                      source, 0.05)
    colormap_vertex_color('../data/field.off',
                          mesh.points(),
                          mesh.face_vertex_indices(),
                          field)

    field = pygeodesic.distance_field(mesh.points(),
                                      mesh.face_vertex_indices(),
                                      source, edge_split=0.05,
                                      max_propagation_distance=1.0)
    field[field>1.e10] = 2.0
    colormap_vertex_color('../data/field_1.off',
                          mesh.points(),
                          mesh.face_vertex_indices(),
                          field)


if __name__ == '__main__':
    # help(pygeodesic)
    test_geodesic()

