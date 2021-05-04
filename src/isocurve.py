"""
=========
ISO-CURVE
=========

Iso-curve extraction from a scalar field defined on a triangular mesh.
A C++ version will be available.

Author: Zishun Liu <liuzishun@gmail.com>
Date: Feb 13, 2021
"""

import numpy as np
import openmesh as om
import networkx as nx
from enum import Enum

__all__ = ['IsoCurve']


class CurveSegmentDir(Enum):
    Flexible = 0
    Determined = 1


class IsoCurve:
    def __init__(self, points, faces, field):
        """
        Iso-curve extraction from a scalar field defined on a triangular mesh.

        Parameters
        ----------
        points : m x 3 array (float)
            mesh vertices
        faces : n x 3 array (int)
            mesh faces
        field : m vector (float)
            scalar field values on mesh vertices
        """

        assert(points.shape[0] == field.size)
        self.mesh = om.TriMesh(points, faces)
        self.field = field.copy()
        self.V = self.mesh.points()
        self.F = self.mesh.face_vertex_indices()

    def extract(self, val, eqlTol=1.e-6):
        """
        Extract the isocurve of value val.

        Parameters
        ----------
        val : float
            extract the iso-curve of value val
        eqlTol : float
             iso-curve passes vertex v if the value on v is within the
             tolerance eqlTol of val
        Returns
        -------
        pts : m x 3 numpy.array (float)
        on_edge : m x 2 numpy.array (int)
        ratio : m numpy.array (float)
        curve_all : list of list
        """

        self.idxOnVerts = np.full(self.mesh.n_vertices(), -1)
        self.idxOnEdges = np.full(self.mesh.n_edges(), -1)

        # init self.pts, self.on_edge, self.ratio
        self._collect_pts(val, eqlTol)

        self.curve_segments = []
        self.curve_segments_dir = []
        for i in range(self.F.shape[0]):
            f = self.F[i]
            cnt = np.count_nonzero(self.idxOnVerts[f] >= 0)
            if cnt == 3:
                print("""[Warning!] Not handled: a whole triangle lays on
                the iso value.""")
            elif cnt == 2:
                self._extract_on_face_through2(i, val)
            elif cnt == 1:
                self._extract_on_face_through1(i, val)
            elif cnt == 0:
                self._extract_on_face_through0(i, val)

        # result topology
        self._construct_topology()

        return self.pts, self.on_edge, self.ratio, self.curve_all

    def _collect_pts(self, val, eqlTol):
        # on verts
        mask = (np.abs(self.field-val) < eqlTol)
        cntPts = np.count_nonzero(mask)
        self.idxOnVerts[mask] = np.arange(cntPts)

        idx = np.where(mask)[0]
        self.pts = [self.V[i] for i in idx]
        self.on_edge = [[i, i] for i in idx]
        self.ratio = [0.0] * cntPts

        # on edges
        for eh in self.mesh.edges():
            heh = self.mesh.halfedge_handle(eh, 0)
            vh0 = self.mesh.from_vertex_handle(heh)
            vh1 = self.mesh.to_vertex_handle(heh)
            if (self.idxOnVerts[vh0.idx()] >= 0 or
                    self.idxOnVerts[vh1.idx()] >= 0):
                continue
            x0 = self.field[vh0.idx()]
            x1 = self.field[vh1.idx()]
            if (x0 < val and x1 > val) or (x0 > val and x1 < val):
                t0 = np.abs(x0-val)
                t1 = np.abs(x1-val)
                t = t0+t1
                t0 /= t
                t1 /= t
                self.pts.append(t1*self.V[vh0.idx()]+t0*self.V[vh1.idx()])
                self.on_edge.append([vh0.idx(), vh1.idx()])
                self.ratio.append(t0)
                self.idxOnEdges[eh.idx()] = cntPts
                cntPts += 1

        self.pts = np.array(self.pts)
        self.on_edge = np.array(self.on_edge)
        self.ratio = np.array(self.ratio)
        return

    def _extract_on_face_through2(self, i, val):
        fh = self.mesh.face_handle(i)
        # find heh
        for h in self.mesh.fh(fh):
            vh0 = self.mesh.from_vertex_handle(h)
            vh1 = self.mesh.to_vertex_handle(h)
            if (self.idxOnVerts[vh0.idx()] >= 0 and
                    self.idxOnVerts[vh1.idx()] >= 0):
                heh = h
                break

        #       vh1
        #        +
        #       /|\
        #      / | \
        # vh2 +--+--+ vh3
        #       vh0

        heh2 = self.mesh.next_halfedge_handle(heh)
        vh2 = self.mesh.to_vertex_handle(heh2)
        if (self.field[vh2.idx()] < val):
            self.curve_segments.append([self.idxOnVerts[vh0.idx()],
                                        self.idxOnVerts[vh1.idx()]])
        else:  # self.field[vh2.idx()] > val. Cannot be == val
            self.curve_segments.append([self.idxOnVerts[vh1.idx()],
                                        self.idxOnVerts[vh0.idx()]])

        heh1 = self.mesh.opposite_halfedge_handle(heh)
        if self.mesh.is_boundary(heh1):
            self.curve_segments_dir.append(CurveSegmentDir.Determined)
            return
        # not boundary
        heh3 = self.mesh.next_halfedge_handle(heh1)
        vh3 = self.mesh.to_vertex_handle(heh3)
        if (self.idxOnVerts[vh3.idx()] or
                (self.field[vh2.idx()] < val and self.field[vh3.idx()] > val) or
                (self.field[vh2.idx()] > val and self.field[vh3.idx()] < val)):
            self.curve_segments_dir.append(CurveSegmentDir.Determined)
            return

        print("[Warning!] Local extrema edge!")
        self.curve_segments_dir.append(CurveSegmentDir.Flexible)
        return

    def _extract_on_face_through1(self, i, val):
        fh = self.mesh.face_handle(i)
        for heh in self.mesh.fh(fh):
            vh = self.mesh.from_vertex_handle(heh)
            if self.idxOnVerts[vh.idx()] < 0:
                continue
            heh1 = self.mesh.next_halfedge_handle(heh)
            eh1 = self.mesh.edge_handle(heh1)
            if self.idxOnEdges[eh1.idx()] < 0:
                # only passes a vertex, does not enter this triangle.
                return
            vh1 = self.mesh.to_vertex_handle(heh)
            if self.field[vh1.idx()] > val:
                self.curve_segments.append([self.idxOnVerts[vh.idx()],
                                            self.idxOnEdges[eh1.idx()]])
                self.curve_segments_dir.append(CurveSegmentDir.Determined)
            else:
                self.curve_segments.append([self.idxOnEdges[eh1.idx()],
                                            self.idxOnVerts[vh.idx()]])
                self.curve_segments_dir.append(CurveSegmentDir.Determined)
            return

    def _extract_on_face_through0(self, i, val):
        fh = self.mesh.face_handle(i)
        for heh in self.mesh.fh(fh):
            eh = self.mesh.edge_handle(heh)
            vh = self.mesh.from_vertex_handle(heh)
            if (self.idxOnEdges[eh.idx()] < 0 or self.field[vh.idx()] > val):
                continue
            heh1 = self.mesh.next_halfedge_handle(heh)
            eh1 = self.mesh.edge_handle(heh1)
            if self.idxOnEdges[eh1.idx()] >= 0:
                self.curve_segments.append([self.idxOnEdges[eh.idx()],
                                            self.idxOnEdges[eh1.idx()]])
                self.curve_segments_dir.append(CurveSegmentDir.Determined)
                return
            # actually here should have an "else"
            heh1 = self.mesh.next_halfedge_handle(heh1)
            eh1 = self.mesh.edge_handle(heh1)
            if self.idxOnEdges[eh1.idx()] >= 0:
                self.curve_segments.append([self.idxOnEdges[eh.idx()],
                                            self.idxOnEdges[eh1.idx()]])
                self.curve_segments_dir.append(CurveSegmentDir.Determined)
                return
            # indeed should have an "else"
            assert False, "Isocurve cuts through one edge. Not possible!"

    def _construct_topology(self):
        num_pts = self.pts.shape[0]

        # construct graph
        G = nx.empty_graph(num_pts)
        for i in range(len(self.curve_segments)):
            if self.curve_segments_dir[i] == CurveSegmentDir.Determined:
                G.add_edge(self.curve_segments[i][0],
                           self.curve_segments[i][1],
                           from_=self.curve_segments[i][0])
            else:
                G.add_edge(self.curve_segments[i][0],
                           self.curve_segments[i][1],
                           from_=-1)

        # check manifoldness
        for i in range(num_pts):
            d = G.degree[i]
            if d == 0:
                print("[Warning!] Isolated point.")
            elif d == 1:
                pass
            elif d == 2:
                neighbors = [j for j in G[i]]
                from0 = G[i][neighbors[0]]['from_']
                from1 = G[i][neighbors[1]]['from_']
                if from0 == from1 or (from0 != i and from1 != i):
                    # both in or both out
                    print("[Warning!] Non-manifold curve.")
            else:
                print("[Warning!] Non-manifold. more than 2 neighbors.")

        self.curve_all = []
        pts_visited = np.zeros(num_pts)
        start = self._find_zero(pts_visited)
        while start >= 0:
            curve = [start]
            pts_visited[start] = 1
            invert = self._extend(curve, G, pts_visited)
            if invert == 1:
                curve = curve[::-1]
            self.curve_all.append(curve)

            start = self._find_zero(pts_visited)

    def _find_zero(self, vec):
        found = np.where(vec == 0)[0]
        if found.size > 0:
            return found[0]
        return -1

    def _extend(self, curve, G, pts_visited):
        # return value:
        # 1: invert 0: flexible -1: not invert.

        neighbors = [i for i in G[curve[0]]]
        if len(neighbors) == 0:
            # curve[0] already marked visited.
            return 0

        # not really front/back
        invert = self._extend_front(curve, neighbors[0], G, pts_visited)
        if len(neighbors) == 2 and curve[0] != curve[-1]:
            # i.e. not a ring
            invert1 = self._extend_back(curve, neighbors[1], G, pts_visited)
            if invert1 + invert == 0 and invert != 0:
                # 1 and -1
                print("[Warning!] Complex direction!")
            invert = max(invert, invert1)

        return invert

    def _extend_front(self, curve, front, G, visited):
        invert = 0

        prev = curve[0]
        now = front
        if visited[now] == 1:
            return invert
        while True:
            visited[now] = 1
            curve.append(now)
            found_next = False
            for n in G[now]:
                if n != prev:
                    prev = now
                    now = n
                    found_next = True
                    if G[prev][now]['from_'] == prev:
                        if (invert == 1):
                            print("[Warning!] Complex direction!")
                        invert = -1  # not invert
                    if G[prev][now]['from_'] == now:
                        if (invert == -1):
                            print("[Warning!] Complex direction!")
                        invert = 1
                    break
            if not found_next:
                break
            if visited[now] == 1:  # ring
                curve.append(now)
                break
        return invert

    def _extend_back(self, curve, back, G, visited):
        # must not be a ring

        invert = 0

        prev = curve[0]
        now = back
        while visited[now] != 1:
            visited[now] = 1
            curve.insert(0, now)
            for n in G[now]:
                if n != prev:
                    prev = now
                    now = n
                    if G[prev][now]['from_'] == prev:
                        if invert == -1:
                            print("[Warning!] Complex direction!")
                        invert = 1
                    if G[prev][now]['from_'] == now:
                        if invert == 1:
                            print("[Warning!] Complex direction!")
                        invert = -1
                    break
        return invert


if __name__ == '__main__':
    pass
