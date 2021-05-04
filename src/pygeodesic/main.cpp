// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.
#include "mmp_wrapper.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

std::tuple<Eigen::MatrixXi, Eigen::VectorXd>
find_path(const Eigen::MatrixXd& V,
          const Eigen::MatrixXi& F,
          int src,
          int dst)
{
    Eigen::MatrixXi path_edge;
    Eigen::VectorXd path_ratio;
    mmp_wrapper wrapper;
    wrapper.find_path(V, F, src, dst, path_edge, path_ratio);
    return std::make_tuple(path_edge, path_ratio);
}

Eigen::VectorXd
distance_field(const Eigen::MatrixXd& V,
               const Eigen::MatrixXi& F,
               const Eigen::VectorXi& source_indices,
               double edge_split = -1.0)
{
    Eigen::VectorXd D;

    mmp_wrapper wrapper;
    wrapper.distance_field(V, F, source_indices, D, edge_split);
    return D;
}

PYBIND11_MODULE(pygeodesic, m) {
    m.doc() = "Python binding of MMP geodesic library";

    // todo: add help for functions

    m.def("find_path", &find_path, "");
    m.def("distance_field", &distance_field, "");
}
