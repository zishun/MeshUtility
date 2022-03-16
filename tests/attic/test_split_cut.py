import openmesh as om
import meshutility as mu
import numpy as np
from numpy.polynomial import Chebyshev as T


def test_cut():
    mesh = om.read_trimesh('../data/sphere_cvt.obj')
    field = mesh.points()[:, 1] - mesh.points()[:, 2]

    isocurve = mu.IsoCurve(mesh.points(), mesh.face_vertex_indices(), field)
    pts, on_edges, ratios, isocurve_indices = isocurve.extract(0.)
    mu.write_obj_lines('../data/curve.obj', pts, isocurve_indices)
    # print(on_edges, ratios)
    mesh, curve_idx = mu.split_mesh(mesh.points(),
                                    mesh.face_vertex_indices(),
                                    on_edges, ratios)
    # print(curve_idx)
    om.write_mesh('../data/mesh_split.obj', mesh)
    curve_idx.append(curve_idx[0])
    curve = [curve_idx[v] for v in isocurve_indices[0]]
    mesh, curve_idx = mu.cut_along_curve(mesh.points(),
                                         mesh.face_vertex_indices(),
                                         curve)
    om.write_mesh('../data/mesh_cut.obj', mesh)

    parts = mu.split_connected_components(mesh.points(),
                                          mesh.face_vertex_indices())
    for i in range(len(parts)):
        om.write_mesh('../data/part_%d.obj' % i, parts[i])


def test_split_rect_0():
    mesh = om.read_trimesh('../data/rect.obj')
    pts = mesh.points()
    field = pts[:, 0] - pts[:, 1]

    isocurve = mu.IsoCurve(pts, mesh.fv_indices(), field)
    pts, on_edges, ratios, isocurve_indices = isocurve.extract(0.)
    mu.write_obj_lines('../data/curve.obj', pts, isocurve_indices)
    # print(on_edges, ratios)
    mesh, curve_idx = mu.split_mesh(mesh.points(),
                                    mesh.face_vertex_indices(),
                                    on_edges, ratios)
    # print(curve_idx)
    om.write_mesh('../data/mesh_split.obj', mesh)
    

def test_split_rect_1():
    # two functions from https://stackoverflow.com/a/9997374
    def ccw(A,B,C):
        return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

    # Return true if line segments AB and CD intersect
    def intersect(A,B,C,D):
        return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

    mesh = om.read_trimesh('../data/rect.obj')
    mesh_pts = mesh.points()
    mesh_edges = mesh.ev_indices()
    
    x = np.linspace(-1, 1, 101)
    y = T.basis(5)(x)
    
    # export Chebyshev polynomial
    pts = np.zeros((x.size, 3))
    pts[:, 0] = x
    pts[:, 1] = y
    mu.write_obj_lines('../data/curve_Chebyshev.obj', pts, [np.arange(pts.shape[0])])

    edges = [np.array([0,0], 'i')]
    ratios = [0.0]
    for i in range(1, x.size-2):
        A = np.array([x[i], y[i]])
        B = np.array([x[i+1], y[i+1]])
        temp_x = []
        temp_e = []
        temp_u = []
        for e in mesh_edges:
            C = mesh_pts[e[0], :2]
            D = mesh_pts[e[1], :2]
            if intersect(A, B, C, D):
                # A = (x1, y1), B = (x2, y2)    
                # C = (x3, y3), D = (x4, y4)
                u = ((A[0]-C[0])*(A[1]-B[1])-(A[1]-C[1])*(A[0]-B[0]))/((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))
                temp_e.append(e)
                temp_u.append(u)
                temp_x.append(C[0]*(1-u)+D[0]*u)
        order = np.argsort(np.array(temp_x))
        edges.extend([temp_e[i] for i in order])
        ratios.extend([temp_u[i] for i in order])
    edges.append(np.array([2,2], 'i'))
    ratios.append(0.0)
    edges = np.array(edges, 'i')
    #print(ratios)
    ratios = np.array(ratios, 'd')
    
    pts0 = mesh_pts[edges[:, 0], :]
    pts1 = mesh_pts[edges[:, 1], :]
    pts = np.multiply(pts0, 1.-ratios[:, np.newaxis]) + \
        np.multiply(pts1, ratios[:, np.newaxis])
    mu.write_obj_lines('../data/curve.obj', pts, [np.arange(pts.shape[0])])
    # print(on_edges, ratios)
    mesh, curve_idx = mu.split_mesh_complete(mesh.points(),
                                    mesh.face_vertex_indices(),
                                    edges, ratios)
    # print(curve_idx)
    om.write_mesh('../data/mesh_split.obj', mesh)
    

if __name__ == '__main__':
    test_cut()
    # test_split_rect_1()
