import openmesh as om

import sys
sys.path.append('../../')
from MeshUtility import IsoCurve, write_obj_lines, split_mesh, \
    cut_along_curve, split_connected_components


def test_cut():
    mesh = om.read_trimesh('../data/sphere_cvt.obj')
    field = mesh.points()[:, 2]

    isocurve = IsoCurve(mesh.points(), mesh.face_vertex_indices(), field)
    pts, on_edges, ratios, isocurve_indices = isocurve.extract(0.)
    write_obj_lines('../data/curve.obj', pts, isocurve_indices)
    # print(on_edges, ratios)
    mesh, curve_idx = split_mesh(mesh.points(),
                                 mesh.face_vertex_indices(),
                                 on_edges, ratios)
    # print(curve_idx)
    om.write_mesh('../data/mesh_split.obj', mesh)
    curve_idx.append(curve_idx[0])
    curve = [curve_idx[v] for v in isocurve_indices[0]]
    mesh, curve_idx = cut_along_curve(mesh.points(),
                                      mesh.face_vertex_indices(),
                                      curve)
    om.write_mesh('../data/mesh_cut.obj', mesh)

    parts = split_connected_components(mesh.points(),
                                       mesh.face_vertex_indices())
    for i in range(len(parts)):
        om.write_mesh('../data/part_%d.obj' % i, parts[i])


if __name__ == '__main__':
    test_cut()
