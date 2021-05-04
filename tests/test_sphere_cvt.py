import openmesh as om

import sys
sys.path.append('../../')
from MeshUtility import sphere_cvt

def test_sphere_cvt():
    verts, faces = sphere_cvt(100, iters=100)
    mesh = om.TriMesh(verts, faces)
    om.write_mesh('../data/sphere_cvt.obj', mesh)


if __name__ == '__main__':
    test_sphere_cvt()
