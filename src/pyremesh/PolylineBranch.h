#ifndef POLYLINE_BRANCH_H
#define POLYLINE_BRANCH_H

#include <vector>
#include <Eigen/Dense>


class PolylineBranch{
/* PolylineBranch is a simple 1D manifold polyline curve (might be a closed curve),
 * It's ends may be *strong* features, usually corners/T-joints etc.
 * Each interior vertex has a unique parameter between 0 and 1.
 * The ends do not have a unique parameter, set as -1.
 */
public:
    PolylineBranch(const std::vector<int>& idx);
    ~PolylineBranch(){}

    void project_to_branch(const Eigen::Vector3d& p, int t0, int t1, int t2,
                           const Eigen::MatrixX3d& V,
                           Eigen::Vector3d& res_p, int& res_t);


private:
    // to remove
    void project_to_branch_test(const Eigen::Vector3d& p, int t0, int t1, int t2,
                           const Eigen::MatrixX3d& V,
                           Eigen::Vector3d& res_p, int& res_t);

    void project_to_line_segment(const Eigen::Vector3d& p,
                                 const Eigen::Vector3d& p0,
                                 const Eigen::Vector3d& p1,
                                 Eigen::Vector3d& res);

public:
    // idx_[0] == idx_[-1] if it is a loop.
    std::vector<int> idx_;
    int N_;
};


#endif  // POLYLINE_BRANCH_H defined
