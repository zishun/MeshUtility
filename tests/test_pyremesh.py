import numpy as np
import openmesh as om
import meshutility as mu


def rect():
    pts = np.array([
        [-1, -1, 0],
        [ 1, -1, 0],
        [ 1,  1, 0],
        [-1,  1, 0],
    ], 'f')
    pts[:, 0] += 0.2
    faces = np.array([
        [0, 1, 2],
        [0, 2, 3]
    ], 'i')

    remesher = mu.pyremesh.remesher()
    remesher.init_mesh(pts, faces)
    remesher.set_feature_vertices_by_angle(20)
    remesher.remesh(0.1, 10)
    verts, faces = remesher.get_mesh()
    mesh = om.TriMesh(verts, faces)
    om.write_mesh(folder+'rect.obj', mesh)


def hemisphere():
    def shape(x, y):
        r1_sqr = np.power(x, 2) + np.power(y, 2)
        z_sqr = r*r - r1_sqr
        z_sqr[z_sqr < 0.0] = 0.0
        return np.sqrt(z_sqr)

    r = 1.0
    n = 100
    pts = np.zeros((n+1, 3))
    angles = np.arange(n)*(2*np.pi/n)
    pts[1:, 0] = r*np.cos(angles)
    pts[1:, 1] = r*np.sin(angles)
    faces = np.zeros((n, 3), 'i')
    faces[:, 1] = np.arange(n, dtype='i')+1
    faces[:, 2] = faces[:, 1]+1
    faces[-1, 2] = 1

    remesher = mu.pyremesh.remesher()
    remesher.init_mesh(pts, faces)
    remesher.set_feature_vertices_by_angle(20)

    targetEdgeLength = 2*np.pi*r/n
    low = (4.0 / 5.0) * targetEdgeLength
    high = (4.0 / 3.0) * targetEdgeLength
    for _ in range(100):
        remesher.split_long_edges(high)
        remesher.collapse_short_edges(low, high)
        remesher.equalize_valences()
        remesher.tangential_relaxation()
        
        # project_to_surface()
        verts, _ = remesher.get_mesh()
        verts[:, 2] = shape(verts[:, 0], verts[:, 1])
        remesher.update_points(verts)

    verts, faces = remesher.get_mesh()
    mesh = om.TriMesh(verts, faces)
    om.write_mesh(folder+'hemisphere.obj', mesh)


def sphere():
    # tetrahedron
    pts = np.array([[1, -1, -1],
        [0, 1, -1],
        [-1, -1, -1],
        [0, 0, 1]])
    faces = np.array([[2, 1, 0],
        [0, 1, 3],
        [1, 2, 3],
        [2, 0, 3]], 'i')

    # remesh to sphere
    remesher = mu.pyremesh.remesher()
    remesher.init_mesh(pts, faces)

    targetEdgeLength = 0.2
    low = (4.0 / 5.0) * targetEdgeLength
    high = (4.0 / 3.0) * targetEdgeLength
    for _ in range(10):
        remesher.split_long_edges(high)
        remesher.collapse_short_edges(low, high)
        remesher.equalize_valences()
        remesher.tangential_relaxation()
        
        # project_to_surface()
        verts, _ = remesher.get_mesh()
        norms = np.linalg.norm(verts, axis=1)
        remesher.update_points(verts / norms[:, np.newaxis])

    verts, faces = remesher.get_mesh()
    mesh = om.TriMesh(verts, faces)
    om.write_mesh(folder+'sphere.obj', mesh)


def sphere_with_feature():
    fn_mesh = folder+'sphere.obj'
    mesh = om.read_trimesh(fn_mesh)

    # split
    field = mesh.points()[:, 0]
    
    isocurve = mu.pyisocurve.isocurve(mesh.points(), 
                                      mesh.face_vertex_indices(), 
                                      field)
    pts, on_edges, ratios, isocurve_indices = isocurve.extract(0.0, 1.e-5)
    mu.write_obj_lines(folder+'curve.obj', 
            pts, isocurve_indices)
    mesh, curve_idx = mu.split_mesh(mesh.points(),
                                 mesh.face_vertex_indices(),
                                 on_edges, ratios)
    # print(curve_idx)
    # om.write_mesh('../data/mesh_split.obj', mesh)
    curve = [curve_idx[v] for v in isocurve_indices[0]]
    feature = np.empty((len(curve)-1, 2), 'i')
    feature[:, 0] = curve[:-1]
    feature[:, 1] = curve[1:]

    remesher = mu.pyremesh.remesher()
    remesher.init_mesh(mesh.points(),
                       mesh.face_vertex_indices())
    remesher.set_features(feature)

    feat_edges = remesher.get_features()
    feat_verts = remesher.get_feature_vertices()
    # write_obj_lines(folder+'feature_detect_e.obj', mesh.points(), feat_edges)
    # write_obj_lines(folder+'feature_detect_v.obj', mesh.points()[feat_verts], [])

    remesher.remesh(0.1, 15)
    verts, faces = remesher.get_mesh()
    mesh = om.TriMesh(verts, faces)
    om.write_mesh(folder+'sphere_split.obj', mesh)

    # feat_edges = remesher.get_features()
    # write_obj_lines(folder+'feature_detect_e.obj', mesh.points(), feat_edges)


if __name__ == '__main__':
    folder = '../data/'
    rect()
    hemisphere()
    sphere()
    sphere_with_feature()
