// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.
#include <iostream>

#include <Eigen/Sparse>

#include "isocurve.h"

isocurve::isocurve(const Eigen::MatrixXd& V,
                   const Eigen::MatrixXi& F,
                   const Eigen::VectorXd& field)
{
    if (V.rows() != field.size())
        throw std::invalid_argument( "V.rows() != field.size()" );
    if (F.cols() != 3)
        throw std::invalid_argument( "Not a triangular mesh" );

    V_ = V;
    F_ = F;
    field_ = field;

    mesh_.clean();
    std::vector<EigenTriMesh::VertexHandle> vhandles(V.rows());
    for (Eigen::Index i = 0; i < V.rows(); ++i)
    {
        //vhandles[i] = mesh_.add_vertex(EigenTriMesh::Point(V(i,0), V(i,1), V(i,2)));
        vhandles[i] = mesh_.add_vertex(V.row(i));
    }

    std::vector<EigenTriMesh::VertexHandle> face_vhandles(3);
    for (Eigen::Index i = 0; i < F.rows(); ++i)
    {
        face_vhandles[0] = vhandles[F(i,0)];
        face_vhandles[1] = vhandles[F(i,1)];
        face_vhandles[2] = vhandles[F(i,2)];
        mesh_.add_face(face_vhandles);
    }
}

isocurve::~isocurve(){}

void isocurve::extract(const double val, const double eqlTol, const bool return_path,
                       Eigen::MatrixXd& _pts,
                       Eigen::MatrixXi& _on_edge,
                       Eigen::VectorXd& _ratio,
                       std::vector<std::deque<int>>& _curve_all) const
{
    std::vector<Eigen::Vector3d> pts;
    std::vector<Eigen::Vector2i> on_edge;
    std::vector<double> ratio;
    std::vector<Eigen::Vector2i> curve_segments;
    std::vector<CurveSegmentDir> curve_segments_dir;

    Eigen::VectorXi idxOnVerts(mesh_.n_vertices());
    Eigen::VectorXi idxOnEdges(mesh_.n_edges());
    idxOnVerts.fill(-1);
    idxOnEdges.fill(-1);

    collect_pts(idxOnVerts, idxOnEdges, pts, on_edge, ratio, val, eqlTol);

    for (Eigen::Index i = 0; i < F_.rows(); ++i) {
        Eigen::VectorXi fv = F_.row(i);
        int cnt = 0;
        for (int j = 0; j < 3; ++j) {
            if (idxOnVerts[fv[j]] >= 0)
                cnt++;
        }
        switch (cnt) {
        case 3:
            std::cout << "[Warning!] Not handled: a whole triangle lays on the iso value." << std::endl;
            break;
        case 2:
            extract_on_face_through2(static_cast<unsigned int>(i), val,
                                     idxOnVerts, idxOnEdges,
                                     curve_segments, curve_segments_dir);
            break;
        case 1:
            extract_on_face_through1(static_cast<unsigned int>(i), val,
                                     idxOnVerts, idxOnEdges,
                                     curve_segments, curve_segments_dir);
            break;
        case 0:
            extract_on_face_through0(static_cast<unsigned int>(i), val,
                                     idxOnVerts, idxOnEdges,
                                     curve_segments, curve_segments_dir);
            break;
        default:
            break;
        }
    }

    // collect results
    _pts.resize(pts.size(), 3);
    for (size_t i = 0; i < pts.size(); ++i) {
        _pts.row(i) = pts[i];
    }
    _on_edge.resize(on_edge.size(), 2);
    for (size_t i = 0; i < on_edge.size(); ++i) {
        _on_edge.row(i) = on_edge[i];
    }
    _ratio = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ratio.data(), ratio.size());

    if (return_path) {
        // collect results: topology of paths
        construct_topology(_curve_all, _pts.rows(), curve_segments, curve_segments_dir);
    } else {
        // if paths are not required, return a collection of edges.
        // the last column is the direction info.
        _curve_all.clear();
        for (size_t i = 0; i < curve_segments.size(); ++i) {
            std::deque<int> curve;
            curve.push_back(curve_segments[i][0]);
            curve.push_back(curve_segments[i][1]);
            switch (curve_segments_dir[i]) {
            case CurveSegmentDir::Determined:
                curve.push_back(1);
                break;
            case CurveSegmentDir::Flexible:
                curve.push_back(2);
                break;
            }
            _curve_all.push_back(curve);
        }
    }
}

void isocurve::collect_pts(Eigen::VectorXi& idxOnVerts,
                           Eigen::VectorXi& idxOnEdges,
                           std::vector<Eigen::Vector3d>& pts,
                           std::vector<Eigen::Vector2i>& on_edge,
                           std::vector<double>& ratio,
                           const double val, const double eqlTol) const {
    int cntPts = 0;

    for (Eigen::Index i = 0; i < V_.rows(); ++i) {
        if (std::abs(field_[i]-val) < eqlTol) {
            idxOnVerts[i] = cntPts;
            cntPts++;
            pts.push_back(V_.row(i));
            on_edge.push_back(Eigen::Vector2i(i, i));
            ratio.push_back(0.0);
        }
    }

    // std::vector<bool> face_visited(mesh_.n_faces(), false);  // actually useless

    for (EigenTriMesh::EdgeHandle eh : mesh_.all_edges()) {
        EigenTriMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
        EigenTriMesh::VertexHandle vh0 = mesh_.from_vertex_handle(heh);
        EigenTriMesh::VertexHandle vh1 = mesh_.to_vertex_handle(heh);
        if (idxOnVerts[vh0.idx()] >= 0 || idxOnVerts[vh1.idx()] >= 0) {
            continue;
        }
        double x0 = field_[vh0.idx()];
        double x1 = field_[vh1.idx()];
        if ( (x0 < val && x1 > val) || (x0 > val && x1 < val) ) {
            double t0 = std::abs(x0-val);
            double t1 = std::abs(x1-val);
            double t = t0+t1;
            t0 /= t;
            t1 /= t;
            pts.push_back(t1*V_.row(vh0.idx())+t0*V_.row(vh1.idx()));
            on_edge.push_back(Eigen::Vector2i(vh0.idx(), vh1.idx()));
            ratio.push_back(t0);
            idxOnEdges[eh.idx()] = cntPts;
            cntPts++;
        }
    }
}

void isocurve::extract_on_face_through2(const unsigned int i,
                                        const double val,
                                        const Eigen::VectorXi& idxOnVerts,
                                        const Eigen::VectorXi& idxOnEdges,
                                        std::vector<Eigen::Vector2i>& curve_segments,
                                        std::vector<CurveSegmentDir>& curve_segments_dir
                                        ) const {
    EigenTriMesh::FaceHandle fh = mesh_.face_handle(i);
    EigenTriMesh::HalfedgeHandle heh;
    EigenTriMesh::VertexHandle vh0;
    EigenTriMesh::VertexHandle vh1;
    // find heh
    for (EigenTriMesh::HalfedgeHandle h : mesh_.fh_range(fh)) {
        vh0 = mesh_.from_vertex_handle(h);
        vh1 = mesh_.to_vertex_handle(h);
        if (idxOnVerts[vh0.idx()] >= 0 && idxOnVerts[vh1.idx()] >= 0) {
            heh = h;
            break;
        }
    }

    //        vh1
    //        /|\
    //       / | \
    //      /  |  \
    // vh2 +---|---+ vh3
    //        vh0

    EigenTriMesh::HalfedgeHandle heh2 = mesh_.next_halfedge_handle(heh);
    EigenTriMesh::VertexHandle vh2 = mesh_.to_vertex_handle(heh2);
    if (field_[vh2.idx()] < val) {
        curve_segments.push_back(Eigen::Vector2i(idxOnVerts[vh0.idx()],
                                 idxOnVerts[vh1.idx()]));
    }else{  // field_[vh2.idx()] > val. Cannot be == val
        curve_segments.push_back(Eigen::Vector2i(idxOnVerts[vh1.idx()],
                                 idxOnVerts[vh0.idx()]));
    }

    EigenTriMesh::HalfedgeHandle heh1 = mesh_.opposite_halfedge_handle(heh);
    if (mesh_.is_boundary(heh1)) {
        curve_segments_dir.push_back(CurveSegmentDir::Determined);
        return;
    }
    // not boundary
    EigenTriMesh::HalfedgeHandle heh3 = mesh_.next_halfedge_handle(heh1);
    EigenTriMesh::VertexHandle vh3 = mesh_.to_vertex_handle(heh3);
    if ( idxOnVerts[vh3.idx()] ||
         (field_[vh2.idx()] < val && field_[vh3.idx()] > val) ||
         (field_[vh2.idx()] > val && field_[vh3.idx()] < val) ) {
        curve_segments_dir.push_back(CurveSegmentDir::Determined);
        return;
    }

    std::cout << "[Warning!] Local extrema edge!" << std::endl;
    curve_segments_dir.push_back(CurveSegmentDir::Flexible);
    return;
}

void isocurve::extract_on_face_through1(const unsigned int i,
                                        const double val,
                                        const Eigen::VectorXi& idxOnVerts,
                                        const Eigen::VectorXi& idxOnEdges,
                                        std::vector<Eigen::Vector2i>& curve_segments,
                                        std::vector<CurveSegmentDir>& curve_segments_dir
                                        ) const {
    auto fh = mesh_.face_handle(i);
    for (EigenTriMesh::HalfedgeHandle heh : mesh_.fh_range(fh)) {
        // auto eh = mesh_.edge_handle(heh);
        auto vh = mesh_.from_vertex_handle(heh);
        if (idxOnVerts[vh.idx()] < 0)
        {
            continue;
        }
        auto heh1 = mesh_.next_halfedge_handle(heh);
        auto eh1 = mesh_.edge_handle(heh1);
        if (idxOnEdges[eh1.idx()] < 0) {
            // only passes a vertex, does not enter this triangle.
            return;
        }
        auto vh1 = mesh_.to_vertex_handle(heh);
        if (field_[vh1.idx()] > val) {
            curve_segments.push_back(Eigen::Vector2i(idxOnVerts[vh.idx()],
                                     idxOnEdges[eh1.idx()]));
            curve_segments_dir.push_back(CurveSegmentDir::Determined);
        } else {
            curve_segments.push_back(Eigen::Vector2i(idxOnEdges[eh1.idx()],
                                     idxOnVerts[vh.idx()]));
            curve_segments_dir.push_back(CurveSegmentDir::Determined);
        }
        return;
    }
}

void isocurve::extract_on_face_through0(const unsigned int i,
                                        const double val,
                                        const Eigen::VectorXi& idxOnVerts,
                                        const Eigen::VectorXi& idxOnEdges,
                                        std::vector<Eigen::Vector2i>& curve_segments,
                                        std::vector<CurveSegmentDir>& curve_segments_dir
                                        ) const {
    auto fh = mesh_.face_handle(i);
    for (EigenTriMesh::HalfedgeHandle heh : mesh_.fh_range(fh)) {
        auto eh = mesh_.edge_handle(heh);
        auto vh = mesh_.from_vertex_handle(heh);
        if (idxOnEdges[eh.idx()] < 0 || field_[vh.idx()] > val) {
            continue;
        }
        auto heh1 = mesh_.next_halfedge_handle(heh);
        auto eh1 = mesh_.edge_handle(heh1);
        if (idxOnEdges[eh1.idx()] >= 0) {
            curve_segments.push_back(Eigen::Vector2i(idxOnEdges[eh.idx()],
                                     idxOnEdges[eh1.idx()]));
            curve_segments_dir.push_back(CurveSegmentDir::Determined);
            return;
        }
        // actually here should have an "else"
        heh1 = mesh_.next_halfedge_handle(heh1);
        eh1 = mesh_.edge_handle(heh1);
        if (idxOnEdges[eh1.idx()] >= 0) {
            curve_segments.push_back(Eigen::Vector2i(idxOnEdges[eh.idx()],
                                     idxOnEdges[eh1.idx()]));
            curve_segments_dir.push_back(CurveSegmentDir::Determined);
            return;
        }
        // indeed should have an "else"
        throw std::invalid_argument( "Isocurve cuts through one edge. Not possible!" );
    }
}

void isocurve::construct_topology(std::vector<std::deque<int>>& curve_all,
                                  const Eigen::Index num_pts,
                                  std::vector<Eigen::Vector2i>& curve_segments,
                                  const std::vector<CurveSegmentDir>& curve_segments_dir) const {

    typedef Eigen::Triplet<int> T;

    // construct graph
    std::vector<T> triplets;
    for (size_t i = 0; i < curve_segments.size(); ++i) {
        int i0 = curve_segments[i][0];
        int i1 = curve_segments[i][1];
        switch (curve_segments_dir[i]) {
        case CurveSegmentDir::Determined:
            triplets.push_back(T(i0, i1, 1));
            triplets.push_back(T(i1, i0, -1));
            break;
        case CurveSegmentDir::Flexible:
            triplets.push_back(T(i0, i1, 2));
            triplets.push_back(T(i1, i0, 2));
            break;
        //default:
        //    break;
        }
    }

    SpMat G(num_pts,num_pts);
    G.setFromTriplets(triplets.begin(), triplets.end());

    // std::vector<int> pts_directed(num_pts, 0);
    // std::vector<int> pts_ends;

    // check manifoldness
    std::vector<int> neighbors;
    std::vector<int> neighbor_tags;
    for (int k=0; k < G.outerSize(); ++k) {
        neighbors.clear();
        neighbor_tags.clear();
        for (SpMat::InnerIterator it(G,k); it; ++it)
        {
            neighbors.push_back(it.col());
            neighbor_tags.push_back(it.value());
//            if (it.value() == 1 || it.value() == -1) {
//                pts_directed[k] = 1;
//            }
        }
        switch (neighbors.size()) {
        case 0:
            std::cout << "[Warning!] Isolated point";
            std::cout << " (#" << k << ")." << std::endl;
            break;
        case 1:
//        if (neighbors.size() == 1) {
//            pts_ends.push_back(k);
//        }
            break;
        case 2:
            if (neighbor_tags[0] == neighbor_tags[1] && std::abs(neighbor_tags[1]) == 1) {
                // both 1 or both -1
                std::cout << "[Warning!] Non-manifold curve";
                std::cout << " (#" << k << ")." << std::endl;
            }
            break;
        default:
            std::cout << "[Warning!] Non-manifold. more than 2 neighbors"; 
            std::cout << " (#" << k << " has " << neighbors.size() << " neighbors)." << std::endl;
            break;
        }
    }
    // sm1.coeffRef(i,j) = v_ij;
    // heads = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(heads_vec.data(), heads_vec.size());

    std::vector<int> pts_visited(num_pts, 0);
    //std::vector<std::deque<int>> curve_all;
    //std::cout << pts_visited[0] << std::endl;
    int start = find_zero(pts_visited);
    //std::cout << "start" << start << std::endl;
    //std::vector<std::deque<int>> curve_all_t;
    while (start >= 0) {
        std::deque<int> curve;
        curve.push_back(start);
        pts_visited[start] = 1;
        int invert = extend(curve, G, pts_visited);
        if (invert == 1)
            std::reverse(curve.begin(), curve.end());
        curve_all.push_back(curve);

        start = find_zero(pts_visited);
    }

//    triplets.clear();
//    for (size_t i = 0; i < curve_all.size(); ++i) {
//        for (size_t j = 0; j < curve_all[i].size()-1; ++j) {
//            triplets.push_back(T(curve_all_t[i][j], curve_all_t[i][j+1], 1));
//            triplets.push_back(T(curve_all_t[i][j], curve_all_t[i][j+1], -1));
//        }
//    }

//    SpMat G1(num_pts,num_pts);
//    G1.setFromTriplets(triplets.begin(), triplets.end());

//    // compare G and G1

}

inline int isocurve::find_zero(const std::vector<int>& vec) const {
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == 0)
            return i;
    }
    return -1;
}

int isocurve::extend(std::deque<int>& curve, const SpMat& G, std::vector<int>& pts_visited) const {
    // return value:
    // 1: invert; 0: flexible; -1: not invert.

    std::vector<int> neighbors;
    //std::vector<int> neighbor_tags;
    for (SpMat::InnerIterator it(G, curve[0]); it; ++it) {
        neighbors.push_back(it.col());
        //neighbor_tags.push_back(it.value());
    }
    if (neighbors.size() == 0) {
        // curve[0] already marked visited.
        return 0;
    }

    // not really front/back
    int invert = extend_front(curve, neighbors[0], G, pts_visited);
    //if (neighbors.size() == 2 && curve[0] != *curve.end()) {
    if (neighbors.size() == 2 && curve[0] != curve.back()) {
        // i.e. not a ring
        int invert1 = extend_back(curve, neighbors[1], G, pts_visited);
        if (invert1 + invert == 0 && invert != 0) {
            // 1 and -1
            printf("[Warning!] Complex direction!");
        }
        invert = std::max(invert, invert1);
    }

    return invert;
}

int isocurve::extend_front(std::deque<int>& curve, int front, const SpMat& G, std::vector<int>& visited) const {
    int invert = 0;

    //int prev = *curve.end();
    int prev = curve[0];
    int now = front;
    if (visited[now] == 1)
        return invert;
    while (true) {
        visited[now] = 1;
        curve.push_back(now);
        bool found_next = false;
        for (SpMat::InnerIterator it(G, now); it; ++it) {
            int n = it.col();
            if (n != prev) {
                prev = now;
                now = n;
                found_next = true;
                if (it.value() == 1) {
                    if (invert == 1)
                        printf("[Warning!] Complex direction!");
                    invert = -1;  // not invert
                }
                if (it.value() == -1) {
                    if (invert == -1)
                        printf("[Warning!] Complex direction!");
                    invert = 1;
                }
                break;
            }
        }
        if (!found_next)
            break;
        if (visited[now] == 1) {  // ring
            curve.push_back(now);
            break;
        }
    }
    return invert;
}

int isocurve::extend_back(std::deque<int>& curve, int back, const SpMat& G, std::vector<int>& visited) const {
    // must not be a ring

    int invert = 0;

    int prev = curve[0];
    int now = back;
    while (visited[now] != 1) {
        visited[now] = 1;
        curve.push_front(now);
        for (SpMat::InnerIterator it(G, now); it; ++it) {
            int n = it.col();
            if (n != prev) {
                prev = now;
                now = n;
                if (it.value() == 1) {
                    if (invert == -1)
                        printf("[Warning!] Complex direction!");
                    invert = 1;
                }
                if (it.value() == -1) {
                    if (invert == 1)
                        printf("[Warning!] Complex direction!");
                    invert = -1;
                }
                break;
            }
        }
    }
    return invert;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, std::vector<std::deque<int>>>
isocurve::extract(const double val, const double eqlTol, const bool return_path) {
    Eigen::MatrixXd pts;
    Eigen::MatrixXi on_edge;
    Eigen::VectorXd ratio;
    std::vector<std::deque<int>> curves_all;
    extract(val, eqlTol, return_path, pts, on_edge, ratio, curves_all);
    // std::cout << curves_all.size() << std::endl;
    auto res = std::make_tuple(pts, on_edge, ratio, curves_all);
    return res;
}
