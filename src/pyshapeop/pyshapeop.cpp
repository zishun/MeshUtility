#include "solverapi.h"

namespace py = pybind11;

PYBIND11_MODULE(pyshapeop, m) {
    m.doc() = "A partial python binding of ShapeOp";

    // todo: add help for functions

    py::class_<solver>(m, "solver")
            .def(py::init())
            .def("setPoints", &solver::setPoints)
            .def("getPoints", &solver::getPoints)
            .def("initialize", &solver::initialize)
            .def("solve", &solver::solve)
            .def("addClosenessConstraint", &solver::addClosenessConstraint)
            .def("addSimilarityConstraint", &solver::addSimilarityConstraint)
            .def("addUniformLaplacianConstraint", &solver::addUniformLaplacianConstraint);

}
