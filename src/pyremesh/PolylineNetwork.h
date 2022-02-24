#ifndef POLYLINE_NETWORK_H
#define POLYLINE_NETWORK_H

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "PolylineBranch.h"


class PolylineNetwork{

    typedef Eigen::SparseMatrix<int, Eigen::RowMajor> SpMat;
    typedef Eigen::Triplet<int> T;

public:
    PolylineNetwork(){}
    ~PolylineNetwork(){}

    void construct_network(const Eigen::MatrixX3d& V,
                           const Eigen::MatrixXi& E,
                           const Eigen::VectorXi& S,
                           Eigen::MatrixXi& Pos);

    // already known that the point p lies on the branch between t0 and t1+1.
    void project_to_network(const Eigen::Vector3d& p,
                            const int branch, int t0, int t1, int t2,
                            const Eigen::MatrixX3d& V, Eigen::Vector3d& res_p, int& res_t);

    int find_branch(int v0, int v1);
    int find_segment(int branch, int v);

    // debug. to be removed.
    void print_network();

private:
    //Eigen::MatrixX3d V_;
    std::vector<std::unique_ptr<PolylineBranch>> branches_;
};


#endif  // POLYLINE_NETWORK_H defined
