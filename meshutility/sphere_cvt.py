"""
==========
SPHERE CVT
==========

Iteratively approximated centroidal Voronoi tessellation (CVT) on the unit
sphere. It samples N points uniformly on the unit sphere.
This is an implmentation of the algorithm described in
https://people.sc.fsu.edu/~jburkardt/m_src/sphere_cvt/sphere_cvt.html

by Zishun Liu <liuzishun@gmail.com>
Feb 13, 2021
"""

import numpy as np
from scipy.spatial import SphericalVoronoi

__all__ = ['sphere_cvt']


def sphere_cvt(n, iters=100):
    """
    Iteratively approximated centroidal Voronoi tessellation (CVT) on the unit
    sphere.

    Parameters
    ----------
    n : int
        # points
    iters : int
        # iterations
    Returns
    -------
    verts : m x 3 array (float)
        Voronoi region centers
    faces : n x 3 array (int)
        Delaunay triangulation of verts.
    """
    scvt = SphereCVT(n, iters)
    verts = scvt.get_points()
    faces = scvt.get_faces()
    return verts, faces


class SphereCVT:
    def __init__(self, n, iters=100):
        self.num = n
        self.points = self.uniform(n)
        for i in range(iters):
            sv = SphericalVoronoi(self.points)
            for j in range(len(sv.regions)):
                region = sv.regions[j]
                p = np.mean(sv.vertices[region, :], axis=0)
                p /= np.linalg.norm(p)
                self.points[j, :] = p[:]

            print('%d/%d\r' % (i+1, iters), end='')
        print()  # close \r
        self.sv = sv

    def uniform(self, num):
        vec = np.random.randn(3, num)
        vec /= np.linalg.norm(vec, axis=0)
        return vec.T

    def get_points(self):
        return self.points

    def get_faces(self):
        num_v = np.max([np.max(x) for x in self.sv.regions])+1
        faces = np.zeros((num_v, 3), dtype=int)
        idx = np.zeros((num_v,), dtype=int)
        for i in range(len(self.sv.regions)):
            for v in self.sv.regions[i]:
                faces[v, idx[v]] = i
                idx[v] += 1
        # re-orient
        for i in range(faces.shape[0]):
            v0 = self.points[faces[i, 0], :]
            e1 = self.points[faces[i, 1], :] - v0
            e2 = self.points[faces[i, 2], :] - v0
            n = np.cross(e1, e2)
            if np.dot(v0, n) < 0:
                faces[i, 1], faces[i, 2] = faces[i, 2], faces[i, 1]
        return faces


if __name__ == '__main__':
    pass
