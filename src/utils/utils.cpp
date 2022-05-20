#include <iostream>
#include <Eigen/Dense>

#include "igl/cut_mesh.h"
#include "igl/barycentric_coordinates.h"
#include "igl/point_mesh_squared_distance.h"

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


Eigen::VectorXd 
get_scalar_field_on_resampled_points(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& field, const Eigen::MatrixXd& P)
{
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(P,V,F,sqrD,I,C);
    
    // Vectorize this?
    int n = P.rows();
    Eigen::MatrixXd bary;
    Eigen::VectorXd result(n);
    for (int i = 0; i < n; ++i)
    {
        igl::barycentric_coordinates(C.row(i), V.row(F(I(i),0)), V.row(F(I(i),1)), V.row(F(I(i),2)), bary);
        result(i) = bary(0)*field(F(I(i),0)) + \
            bary(1)*field(F(I(i),1)) + \
            bary(2)*field(F(I(i),2));
    }
    
    return result;
}


PYBIND11_MODULE(utils, m) {
    m.doc() = "";

    // todo: add help for functions
    // keyword / default arguments: https://pybind11.readthedocs.io/en/stable/basics.html#keyword-arguments

    m.def("igl_cut_mesh", &igl_cut_mesh, "");
    m.def("get_scalar_field_on_resampled_points", &get_scalar_field_on_resampled_points, "");
}
