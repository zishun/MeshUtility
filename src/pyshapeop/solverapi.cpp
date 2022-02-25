#include "solverapi.h"

solver::solver()
{
    s_ = ShapeOp::Solver();
}

solver::~solver(){}

void solver::setPoints(const Eigen::MatrixXd& points)
{
    ShapeOp::Matrix3X mat = points.transpose();
    
    s_.setPoints(mat);
}

void solver::addClosenessConstraint(int i, Eigen::VectorXd& target, double weight)
{
    std::vector<int> id_vector;
    id_vector.push_back(i);
    auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, weight, s_.getPoints());
    ShapeOp::Vector3 position = target;
    c->setPosition(position);
    s_.addConstraint(c);
}

void solver::addSimilarityConstraint(std::vector<int>& face, 
                                     const Eigen::MatrixXd& shape,
                                     double weight,
                                     bool scaling, bool rotate, bool flip)
{
    auto c = std::make_shared<ShapeOp::SimilarityConstraint>(
                        face, weight, s_.getPoints(), scaling, rotate, flip);
    s_.addConstraint(c);

    std::vector<ShapeOp::Matrix3X> shapes;
    ShapeOp::Matrix3X shape0 = shape.transpose();
    
    shapes.push_back(shape0);
    c->setShapes(shapes);
}

void solver::addUniformLaplacianConstraint(std::vector<int>& idI,
                                     double weight, bool displacement_lap)
{
    auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(
                        idI, weight, s_.getPoints(), displacement_lap);
    s_.addConstraint(c);
}


void solver::initialize()
{
    s_.initialize();
}

void solver::solve(unsigned int num_inters)
{
    s_.solve(num_inters);
}

Eigen::MatrixXd solver::getPoints()
{
    ShapeOp::Matrix3X mat = s_.getPoints();
    Eigen::MatrixXd pts = mat.transpose();
    return pts;
}
