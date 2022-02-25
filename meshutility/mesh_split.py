"""
==========
SPLIT MESH
==========

Split mesh.

by Zishun Liu <liuzishun@gmail.com>
Feb 13, 2021
"""

import numpy as np
import openmesh as om

__all__ = ['split_mesh', 'split_mesh_complete']


def split_mesh(V, F, edges, ratios):
    """
    Split mesh.

    Parameters
    ----------
    V : numpy.array
    F : numpy.array
    edges : numpy.array
    ratios : numpy.array
        p0 * (1-r) + p1 * r
    Returns
    -------
    mesh : openmesh.TriMesh
    index : list (int)
    """

    mesh = om.TriMesh(V, F)
    assert mesh.n_vertices() == V.shape[0], "Invalid input mesh"

    points = mesh.points()
    pts0 = points[edges[:, 0], :]
    pts1 = points[edges[:, 1], :]
    curve_pts = np.multiply(pts0, 1.-ratios[:, np.newaxis]) + \
        np.multiply(pts1, ratios[:, np.newaxis])

    new_index = []
    for i in range(edges.shape[0]):
        if edges[i, 0] == edges[i, 1]:
            new_index.append(edges[i, 0])
        else:
            p = curve_pts[i, :]
            heh = mesh.find_halfedge(mesh.vertex_handle(edges[i, 0]),
                                     mesh.vertex_handle(edges[i, 1]))
            eh = mesh.edge_handle(heh)
            vh = mesh.add_vertex(p)
            mesh.split_edge(eh, vh)
            new_index.append(vh.idx())

    return mesh, new_index


def split_mesh_complete(V, F, edges, ratios):
    """
    Split mesh. This is the complete version. Maybe overshoot.

    Parameters
    ----------
    V : numpy.array
    F : numpy.array
    edges : numpy.array, given in the order of the curve
    ratios : numpy.array, corresponding to rows in edges.
        p0 * (1-r) + p1 * r
    Returns
    -------
    mesh : openmesh.TriMesh
    index : list (int)
    """

    # unify the edges/ratios, make sure e[0] <= e[1]
    E = edges.copy()
    R = ratios.copy()
    M = om.PolyMesh(V, F)
    
    M.edge_property('split')
    for eh in M.edges():
        M.set_edge_property('split', eh, [])
    for i in range(R.size):
        if E[i, 0] != E[i, 1]:
            if E[i, 1] < E[i, 0]:
                E[i, 1], E[i, 0] = E[i, 0], E[i, 1]
                R[i] = 1.0 - R[i]
            vh0 = M.vertex_handle(E[i, 0])
            vh1 = M.vertex_handle(E[i, 1])
            heh = M.find_halfedge(vh0, vh1)
            assert heh.is_valid(), 'invalid input edges (row %d)'%i
            eh = M.edge_handle(heh)
            L = M.edge_property('split', eh)
            if len(L) > 0:
                # print('split the same edge several times!')
                L.append(i)
                M.set_edge_property('split', eh, L)
            else:
                M.set_edge_property('split', eh, [i])

    pts0 = V[E[:, 0], :]
    pts1 = V[E[:, 1], :]
    curve_pts = np.multiply(pts0, 1.-R[:, np.newaxis]) + \
        np.multiply(pts1, R[:, np.newaxis])

    new_index = np.empty((R.size,), 'i')
    for i in range(E.shape[0]):
        if E[i, 0] == E[i, 1]:
            new_index[i] = E[i, 0]

    for eh in M.edges():
        L = M.edge_property('split', eh)
        if L is None or len(L) == 0:
            continue
        elif len(L) == 1:
            i = L[0]
            p = curve_pts[i, :]
            heh = M.find_halfedge(M.vertex_handle(E[i, 0]),
                                  M.vertex_handle(E[i, 1]))
            eh_ = M.edge_handle(heh)
            vh = M.add_vertex(p)
            M.split_edge(eh_, vh)
            new_index[i] = vh.idx()
        else:
            heh = M.halfedge_handle(eh, 0)
            vh0 = M.from_vertex_handle(heh)
            vh1 = M.to_vertex_handle(heh)
            if vh0.idx() > vh1.idx():
                vh0, vh1 = vh1, vh0
            L.sort(key=lambda x:R[x])
            for i in L:
                p = curve_pts[i, :]
                vh = M.add_vertex(p)
                heh = M.find_halfedge(vh0, vh1)
                assert(heh.is_valid())
                eh_ = M.edge_handle(heh)
                M.split_edge(eh_, vh)
                vh0 = vh
                new_index[i] = vh.idx()

    def find_shared_face(v0, v1):
        f0 = [fh.idx() for fh in M.vf(v0)]
        f1 = [fh.idx() for fh in M.vf(v1)]
        f = list(set(f0).intersection(f1))
        if len(f) != 1:
            print(len(f))
            assert(False)
        return M.face_handle(f[0])

    for i in range(R.size-1):
        v0 = M.vertex_handle(new_index[i])
        v1 = M.vertex_handle(new_index[i+1])
        heh = M.find_halfedge(v0, v1)
        if not heh.is_valid():
            # find the face containing both vertices
            # split the face by adding edge (v0, v1)
            fh = find_shared_face(v0, v1)
            vhs = [vh for vh in M.fv(fh)]
            vs = [vh.idx() for vh in vhs]
            i0 = vs.index(v0.idx())
            i1 = vs.index(v1.idx())
            if i0 > i1:
                i0, i1 = i1, i0
            f0 = vhs[i0:i1+1]
            f1 = vhs[i1:] + vhs[:i0+1]
            M.delete_face(fh, delete_isolated_vertices=False)
            M.add_face(f0)
            M.add_face(f1)
    M.garbage_collection()
    #return M, new_index  # return polymesh

    def triangulate(M, fh):
        # using deque is an overshoot.
        vhs = [vh for vh in M.fv(fh)]
        n = len(vhs)
        if n == 3:
            return
        pts = np.empty((n, 3))
        for i in range(n):
            pts[i] = M.point(vhs[i])
        left = [i-1 for i in range(n)]
        right = [(i+1)%n for i in range(n)]
        edge_vecs = pts[right] - pts
        norm = np.linalg.norm(edge_vecs, axis=1)
        edge_vecs_norms = edge_vecs / norm[:, np.newaxis]
        minus_cosines = np.sum(edge_vecs_norms[left]*edge_vecs_norms, axis=1)
        if n == 4:
            # handle the simple case seperately.
            k = np.argmax(minus_cosines)
            f0 = [vhs[(k+i)%4] for i in range(3)]
            f1 = [vhs[(k+2+i)%4] for i in range(3)]
            M.delete_face(fh, delete_isolated_vertices=False)
            M.add_face(f0)
            M.add_face(f1)
        else:
            # find three smallest angles
            idx = sorted(np.argpartition(minus_cosines, 3)[:3])
            e0 = [i for i in range(idx[0], idx[1]+1)]
            e1 = [i for i in range(idx[1], idx[2]+1)]
            e2 = [i for i in range(idx[2]-n, idx[0]+1)]
            # rotate to make sure e0 is splitted.
            if len(e0) == 2:
                if len(e1) > 2:
                    e0, e1, e2 = e1, e2, e0
                elif len(e2) > 2:
                    e0, e1, e2 = e2, e0, e1
            M.delete_face(fh, delete_isolated_vertices=False)
            for i in range(1, len(e0)-2):
                M.add_face([vhs[i] for i in [e0[i], e0[i+1], e1[-1]]])
            for i in range(len(e1)-1):
                M.add_face([vhs[i] for i in [e1[i], e1[i+1], e0[-2]]])
            for i in range(len(e2)-1):
                M.add_face([vhs[i] for i in [e2[i], e2[i+1], e0[1]]])

    # May generate degenerated triangles if using the build-in triangulation by converting PolyMesh to TriMesh
    # Tiranulate explicitly!
    for fh in M.faces():
        triangulate(M, fh)
    M.garbage_collection()

    mesh = om.TriMesh(M.points(), M.face_vertex_indices())
    return mesh, new_index


if __name__ == '__main__':
    pass
