"""
========
CUT MESH
========

Cut mesh along curve.

by Zishun Liu <liuzishun@gmail.com>
Feb 13, 2021
"""

import numpy as np
import openmesh as om
import networkx as nx

from .mesh2graph import ff_graph

__all__ = ['cut_along_curve']


def cut_along_curve(V, F, curve_index):
    """
    Cut mesh along curve.

    Parameters
    ----------
    V : numpy.array
    F : numpy.array
    curve_index : list
    Returns
    -------
    mesh : openmesh.TriMesh
    index : list (int)
    """

    mesh = om.TriMesh(V, F)
    assert(len(curve_index) > 1)

    be_loop = (curve_index[0] == curve_index[-1])
    if be_loop:
        cutter = MeshCutterLoop()
        return cutter.cut(mesh, curve_index)
    else:
        cutter = MeshCutter()
        return cutter.cut(mesh, curve_index)


class MeshCutterLoop:
    def __init__(self):
        pass

    def cut(self, mesh, curve_index):
        # based on the connected components of the ff graph
        # the method in MeshCutter is also plausible.

        self.mesh = mesh
        self.curve = curve_index

        self.add_vertices()
        self.construct_graph()
        added_verts = self.update_faces()
        return self.mesh, added_verts

    def add_vertices(self):
        num_v = self.mesh.points().shape[0]
        self.new_idx = -np.ones((num_v,), dtype=np.int32)
        assert np.unique(self.curve).size == len(self.curve)-1, """self-
        intersected loop?"""

        for i in range(len(self.curve)-1):
            vh = self.mesh.vertex_handle(self.curve[i])
            self.mesh.add_vertex(self.mesh.point(vh))
            self.new_idx[self.curve[i]] = num_v+i

    def construct_graph(self):
        self.G = ff_graph(self.mesh)

        self.seed_face = np.zeros((self.mesh.n_faces()), dtype=bool)
        for i in range(len(self.curve)-1):
            vh0 = self.mesh.vertex_handle(self.curve[i])
            vh1 = self.mesh.vertex_handle(self.curve[i+1])
            heh0 = self.mesh.find_halfedge(vh0, vh1)
            assert(heh0.is_valid())
            heh1 = self.mesh.opposite_halfedge_handle(heh0)
            if not (self.mesh.is_boundary(heh0) or
                    self.mesh.is_boundary(heh1)):
                fh0 = self.mesh.face_handle(heh0)
                fh1 = self.mesh.face_handle(heh1)
                self.G.remove_edge(fh0.idx(), fh1.idx())
                self.seed_face[fh0.idx()] = True

    def update_faces(self):
        faces_to_delete = []
        faces_to_add = []
        for cc in nx.connected_components(self.G):
            cc = [x for x in cc]
            if not np.any(self.seed_face[cc]):
                continue
            for i in cc:
                fh = self.mesh.face_handle(i)
                fv = [v.idx() for v in self.mesh.fv(fh)]
                if np.any(self.new_idx[fv] >= 0):
                    faces_to_delete.append(fh)
                    face = []
                    for vid in fv:
                        if self.new_idx[vid] >= 0:
                            vid = self.new_idx[vid]
                        face.append(self.mesh.vertex_handle(vid))
                    faces_to_add.append(face)

        [self.mesh.delete_face(fh) for fh in faces_to_delete]
        self.mesh.garbage_collection()
        for f in faces_to_add:
            self.mesh.add_face(f)

        added_verts = []
        for x in self.curve:
            added_verts.append(self.new_idx[x])
        return added_verts


class MeshCutter:
    def __init__(self,):
        pass

    def cut(self, mesh, curve_index):
        self.mesh = mesh
        self.curve = curve_index

        self.add_vertices()
        self.mark_faces_to_update()
        added_verts = self.update_faces()
        return self.mesh, added_verts

    def add_vertices(self):
        num_v = self.mesh.points().shape[0]
        self.new_idx = -np.ones((num_v,), dtype=np.int32)

        vh = self.mesh.vertex_handle(self.curve[0])
        self.start_on_boundary = self.mesh.is_boundary(vh)
        vh = self.mesh.vertex_handle(self.curve[-1])
        self.end_on_boundary = self.mesh.is_boundary(vh)
        start = 0
        end = len(self.curve)
        if not self.start_on_boundary:
            # start in the interior, do not need to duplicate the first point
            start += 1
        if not self.end_on_boundary:
            # end in the interior, do not need to duplicate the last point
            end -= 1
        for i in range(start, end):
            vh = self.mesh.vertex_handle(self.curve[i])
            self.mesh.add_vertex(self.mesh.point(vh))
            self.new_idx[self.curve[i]] = num_v+(i-start)
        for i in range(1, len(self.curve)-1):
            vh = self.mesh.vertex_handle(self.curve[i])
            assert (not self.mesh.is_boundary(vh)), """interior vertice %d of
            the curve cannot be on the mesh boundary. Consider split the
            curve.""" % self.curve[i]

    def mark_faces_to_update(self):
        self.faces_to_delete = np.zeros((self.mesh.n_faces()), dtype=bool)
        if self.start_on_boundary:
            self.mark_faces_to_update_start_on_boundary()

        for i in range(1, len(self.curve)-1):
            # sector between 2 halfedges
            # i-1 -> i -> (i+1)
            het0 = (self.curve[i-1], self.curve[i])
            het1 = (self.curve[i], self.curve[i+1])
            self.mark_left_faces_in_sector(het0, het1)

        if self.end_on_boundary:
            self.mark_faces_to_update_end_on_boundary()

    def mark_faces_to_update_start_on_boundary(self):
        heh0 = self.mesh.find_halfedge(self.mesh.vertex_handle(self.curve[0]),
                                       self.mesh.vertex_handle(self.curve[1]))
        assert(heh0.is_valid())
        assert(not self.mesh.is_boundary(heh0))
        fh = self.mesh.face_handle(heh0)
        self.faces_to_delete[fh.idx()] = True
        pheh = self.mesh.prev_halfedge_handle(heh0)
        while True:
            oheh = self.mesh.opposite_halfedge_handle(pheh)
            if self.mesh.is_boundary(oheh):
                break
            fh = self.mesh.face_handle(oheh)
            self.faces_to_delete[fh.idx()] = True
            pheh = self.mesh.prev_halfedge_handle(oheh)

    def mark_faces_to_update_end_on_boundary(self):
        heh0 = self.mesh.find_halfedge(self.mesh.vertex_handle(self.curve[-2]),
                                       self.mesh.vertex_handle(self.curve[-1]))
        assert(heh0.is_valid())
        assert(not self.mesh.is_boundary(heh0))
        fh = self.mesh.face_handle(heh0)
        self.faces_to_delete[fh.idx()] = True
        nheh = self.mesh.next_halfedge_handle(heh0)
        while True:
            oheh = self.mesh.opposite_halfedge_handle(nheh)
            if self.mesh.is_boundary(oheh):
                break
            fh = self.mesh.face_handle(oheh)
            self.faces_to_delete[fh.idx()] = True
            nheh = self.mesh.next_halfedge_handle(oheh)

    def update_faces(self,):
        faces_to_add = []
        for i in range(self.faces_to_delete.shape[0]):
            if not self.faces_to_delete[i]:
                continue
            face = []
            fh = self.mesh.face_handle(i)
            for v in self.mesh.fv(fh):
                vid = v.idx()
                if self.new_idx[vid] >= 0:
                    vid = self.new_idx[vid]
                face.append(self.mesh.vertex_handle(vid))
            faces_to_add.append(face)

        for i in range(self.faces_to_delete.shape[0]):
            if not self.faces_to_delete[i]:
                continue
            fh = self.mesh.face_handle(i)
            self.mesh.delete_face(fh)
        self.mesh.garbage_collection()
        for f in faces_to_add:
            self.mesh.add_face(f)

        added_verts = []
        for x in self.curve:
            added_verts.append(self.new_idx[x])
        return added_verts

    # find faces in sector between two halfedges, mark them as to-be-deleted
    def mark_left_faces_in_sector(self, het0, het1):
        heh0 = self.mesh.find_halfedge(self.mesh.vertex_handle(het0[0]),
                                       self.mesh.vertex_handle(het0[1]))
        assert(heh0.is_valid())
        assert(not self.mesh.is_boundary(heh0))
        fh = self.mesh.face_handle(heh0)
        self.faces_to_delete[fh.idx()] = True
        nheh = self.mesh.next_halfedge_handle(heh0)
        while True:
            oheh = self.mesh.opposite_halfedge_handle(nheh)
            nheh = self.mesh.next_halfedge_handle(oheh)
            fh = self.mesh.face_handle(nheh)
            self.faces_to_delete[fh.idx()] = True
            if self.mesh.to_vertex_handle(nheh).idx() == het1[1]:
                break


if __name__ == '__main__':
    pass
