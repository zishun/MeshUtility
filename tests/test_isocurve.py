import openmesh as om

import sys
sys.path.append('../../')
from MeshUtility import IsoCurve, write_obj_lines


def test_isocurve0():
    # run python sphere_cvt.py to generate the mesh
    mesh = om.read_trimesh('../data/sphere_cvt.obj')
    field = mesh.points()[:, 1]

    isocurve = IsoCurve(mesh.points(), mesh.face_vertex_indices(), field)
    pts, _, _, isocurve_indices = isocurve.extract(0.5)
    write_obj_lines('../data/curve.obj', pts, isocurve_indices)


if __name__ == '__main__':
    test_isocurve0()
