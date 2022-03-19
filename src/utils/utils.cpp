#include <iostream>
#include <Eigen/Dense>

#include "igl/cut_mesh.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/eigen.h"

namespace py = pybind11;


std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>
igl_cut_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& cuts)
{
    Eigen::MatrixXd Vn;
    Eigen::MatrixXi Fn;
    igl::cut_mesh(V, F, cuts, Vn, Fn);
        
    auto result = std::make_tuple(Vn, Fn);
    return result;
}

PYBIND11_MODULE(utils, m) {
    m.doc() = "";

    // todo: add help for functions
    // keyword / default arguments: https://pybind11.readthedocs.io/en/stable/basics.html#keyword-arguments

    m.def("igl_cut_mesh", &igl_cut_mesh, "");
}
