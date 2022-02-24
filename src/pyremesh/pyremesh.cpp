#include <iostream>
#include <Eigen/Dense>

#include "IncrementalRemesher.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;


std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>
remesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double targetEdgeLength)
{
    IncrementalRemesher remesher;

    remesher.init_mesh(V, F);
    remesher.collect_features();
    remesher.remesh(targetEdgeLength, 10);

    Eigen::MatrixXd V1;
    Eigen::MatrixXi F1;
    remesher.get_mesh(V1, F1);

    auto result = std::make_tuple(V1, F1);
    return result;
}

PYBIND11_MODULE(pyremesh, m) {
    m.doc() = "";

    // todo: add help for functions
    // keyword / default arguments: https://pybind11.readthedocs.io/en/stable/basics.html#keyword-arguments

    m.def("remesh", &remesh, "");

    py::class_<IncrementalRemesher>(m, "remesher")
            .def(py::init())
            .def("init_mesh", &IncrementalRemesher::init_mesh)
            .def("remesh", &IncrementalRemesher::remesh)
            .def("get_mesh", py::overload_cast<>(&IncrementalRemesher::get_mesh))
            .def("update_points", &IncrementalRemesher::update_points)
            .def("set_features", &IncrementalRemesher::set_features)
            .def("set_feature_vertices", &IncrementalRemesher::set_feature_vertices)
            .def("set_features_by_dihedral_angle", &IncrementalRemesher::set_features_by_dihedral_angle)
            .def("set_feature_vertices_by_angle", &IncrementalRemesher::set_feature_vertices_by_angle)
            .def("get_features", py::overload_cast<>(&IncrementalRemesher::get_features))
            .def("get_feature_vertices", py::overload_cast<>(&IncrementalRemesher::get_feature_vertices))
            .def("collect_features", &IncrementalRemesher::collect_features)
            .def("split_long_edges", &IncrementalRemesher::split_long_edges)
            .def("collapse_short_edges", &IncrementalRemesher::collapse_short_edges)
            .def("equalize_valences", &IncrementalRemesher::equalize_valences)
            .def("tangential_relaxation", &IncrementalRemesher::tangential_relaxation)
            .def("project_to_surface", &IncrementalRemesher::project_to_surface);
}
