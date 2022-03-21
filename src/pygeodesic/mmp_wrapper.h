// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.

#ifndef MESHUTILITY_MMP_WRAPPER_H
#define MESHUTILITY_MMP_WRAPPER_H
// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

#include <iostream>

#include <Eigen/Core>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  OMesh;

#include "geodesic_algorithm_exact.h"


class mmp_wrapper
{
public:
    mmp_wrapper();
    ~mmp_wrapper();

    void find_path(const Eigen::MatrixXd& V,
                   const Eigen::MatrixXi& F,
                   int src,
                   int dst,
                   Eigen::MatrixXi& path_edge,
                   Eigen::VectorXd& path_ratio);
    void distance_field(const Eigen::MatrixXd& V,
                        const Eigen::MatrixXi& F,
                        const Eigen::VectorXi& source_indices,
                        Eigen::VectorXd& D,
                        double edge_split,
                        double max_propagation_distance);

private:
    bool VF_eigenmat_to_stlvector(const Eigen::MatrixXd& V,
                                  const Eigen::MatrixXi& F,
                                  std::vector<double>& points,
                                  std::vector<unsigned>& faces);

    double point_ratio_on_edge(geodesic::Point3D p0, geodesic::Point3D p1,
                               geodesic::Point3D p);

    void find_path(const std::vector<double>& points,
                   const std::vector<unsigned>& faces,
                   int src, int dst,
                   std::vector<geodesic::SurfacePoint>& path,
                   Eigen::MatrixXi& path_edge,
                   Eigen::VectorXd& path_ratio);

    double construct_mesh_for_edge_based_geodesic(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            const Eigen::VectorXi& source_indices,
            std::vector<double>& points,
            std::vector<unsigned>& faces,
            Eigen::MatrixXi& inserted_vertex_info);

    void edge_sourced_geodesic(std::vector<double>& points,
                               const std::vector<unsigned>& faces,
                               int num_dst,
                               int segments,
                               Eigen::MatrixXi& inserted_vertex_info,
                               Eigen::VectorXd& D,
                               double max_propagation_distance);

    void vertex_sourced_geodesic(const std::vector<double>& points,
                                 const std::vector<unsigned>& faces,
                                 const Eigen::VectorXi& source_indices,
                                 Eigen::VectorXd& D,
                                 double max_propagation_distance);
};
#endif  // MESHUTILITY_MMP_WRAPPER_H