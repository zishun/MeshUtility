// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.
#include <vector>
#include <array>
#include <iostream>

#include "isocurve.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;


PYBIND11_MODULE(pyisocurve, m) {
    m.doc() = "isocurves";

    // TODO: add help for functions

    py::class_<isocurve>(m, "isocurve")
                .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXi&, const Eigen::VectorXd&>())
                .def("extract", py::overload_cast<double, double, bool>(&isocurve::extract),
                "",
                py::arg("val"), py::arg("eqlTol"), py::arg("return_path")=true);
}
