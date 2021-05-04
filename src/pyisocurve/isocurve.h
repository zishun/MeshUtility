// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.
#ifndef ANAGEOM_ISOCURVE_H
#define ANAGEOM_ISOCURVE_H

#include <vector>
#include <deque>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include <Eigen/Sparse>

class isocurve {

    struct EigenTraits : OpenMesh::DefaultTraits {
        using Point = Eigen::Vector3d;
        using Normal = Eigen::Vector3d;
        using TexCoord2D = Eigen::Vector2d;
    };
    using EigenTriMesh = OpenMesh::TriMesh_ArrayKernelT<EigenTraits>;
    using SpMat = Eigen::SparseMatrix<int, Eigen::RowMajor>;

    enum CurveSegmentDir {
        Determined = 0,
        Flexible = 1
    };

public:
    isocurve(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& field);
    ~isocurve();
    void extract(const double val, const double eqlTol,
                 Eigen::MatrixXd& _pts,
                 Eigen::MatrixXi& _on_edge,
                 Eigen::VectorXd& _ratio,
                 std::vector<std::deque<int>>& _curve_all) const;
                 // Eigen::VectorXi& _nxt,
                 // Eigen::VectorXi& _heads) const;
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, std::vector<std::deque<int>>>
    extract(const double val, const double eqlTol);

private:
    void collect_pts(Eigen::VectorXi& idxOnVerts,
                     Eigen::VectorXi& idxOnEdges,
                     std::vector<Eigen::Vector3d>& pts,
                     std::vector<Eigen::Vector2i>& on_edge,
                     std::vector<double>& ratio,
                     const double val, const double eqlTol) const;
    void extract_on_face_through2(const unsigned int i,
                                  const double val,
                                  const Eigen::VectorXi& idxOnVerts,
                                  const Eigen::VectorXi& idxOnEdges,
                                  std::vector<Eigen::Vector2i>& curve_segments,
                                  std::vector<CurveSegmentDir>& curve_segments_dir
                                  ) const;
    void extract_on_face_through1(const unsigned int i,
                                  const double val,
                                  const Eigen::VectorXi& idxOnVerts,
                                  const Eigen::VectorXi& idxOnEdges,
                                  std::vector<Eigen::Vector2i>& curve_segments,
                                  std::vector<CurveSegmentDir>& curve_segments_dir
                                  ) const;
    void extract_on_face_through0(const unsigned int i,
                                  const double val,
                                  const Eigen::VectorXi& idxOnVerts,
                                  const Eigen::VectorXi& idxOnEdges,
                                  std::vector<Eigen::Vector2i>& curve_segments,
                                  std::vector<CurveSegmentDir>& curve_segments_dir
                                  ) const;

    void construct_topology(// Eigen::VectorXi& nxt,
                            // Eigen::VectorXi& heads,
                            std::vector<std::deque<int>>& curve_all,
                            const Eigen::Index num_pts,
                            std::vector<Eigen::Vector2i>& curve_segments,
                            const std::vector<CurveSegmentDir>& curve_segments_dir
                            ) const;
    inline int find_zero(const std::vector<int>& vec) const;
    int extend(std::deque<int>& curve, const SpMat& G, std::vector<int>& pts_visited) const;
    int extend_front(std::deque<int>& curve, int front, const SpMat& G, std::vector<int>& pts_visited) const;
    int extend_back(std::deque<int>& curve, int back, const SpMat& G, std::vector<int>& pts_visited) const;

private:
    EigenTriMesh mesh_;
    Eigen::VectorXd field_;
    Eigen::MatrixXd V_;
    Eigen::MatrixXi F_;
};

#endif  // ANAGEOM_ISOCURVE_H
