#include "Solver.h"
#include "Constraint.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class solver
// todo: how about inherit ShapeOp::Solver?
{
public:
    solver();
    ~solver();
    void setPoints(const Eigen::MatrixXd& points);
    void addClosenessConstraint(int i, Eigen::VectorXd& target,
                                double weight);
    void addSimilarityConstraint(std::vector<int>& face, 
                                 const Eigen::MatrixXd& shape,
                                 double weight,
                                 bool scaling, bool rotate, bool flip);
    void addUniformLaplacianConstraint(std::vector<int>& idI,
                                     double weight, bool displacement_lap);
    void initialize();
    void solve(unsigned int num_inters);
    Eigen::MatrixXd getPoints();

private:
    ShapeOp::Solver s_;
};
