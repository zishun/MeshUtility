#include "PolylineBranch.h"

#include <limits>

PolylineBranch::PolylineBranch(const std::vector<int>& idx){
    N_ = idx.size();
    idx_.assign(idx.begin(), idx.end());
}

void PolylineBranch::project_to_branch(const Eigen::Vector3d& p,
                                       int t0, int t1, int t2,
                                       const Eigen::MatrixX3d& V,
                                       Eigen::Vector3d& res_p, int& res_t) {
//    project_to_branch_test(p, t0, t1, t2, V,
//                           res_p, res_t);
//    return;

    int lower = 0;
    int upper = N_-2;

    int N = upper - lower + 1; // t0~t0+1, ..., t1~t1+1
    Eigen::MatrixX3d projections(N, 3);
    Eigen::VectorXd dist(N);
    Eigen::Vector3d proj;
    for (int i = 0; i < N; ++i) {
        project_to_line_segment(p, V.row(idx_[lower+i]), V.row(idx_[lower+i+1]), proj);
        projections.row(i) = proj;
        dist(i) = (proj-p).norm();
    }

    double min_dist = std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i) {
        if (dist(i) < min_dist) {
            min_dist = dist(i);
            res_t = i;
        }
    }
    res_p = projections.row(res_t);
    res_t += lower;

    // checking
//    Eigen::Vector3d res_p_check;
//    int res_t_check;
//    project_to_branch_test(p, t0, t1, t2, V,
//                           res_p_check, res_t_check);
//    if (res_t_check != res_t) {
//        printf("Wrong! got %d, should be %d\n", res_t_check, res_t);
//        printf("v %f %f %f\n", p[0], p[1], p[2]);
//        printf("v %f %f %f\n", V(idx_[t0], 0), V(idx_[t0], 1), V(idx_[t0], 2));
//        printf("v %f %f %f\n", V(idx_[t1], 0), V(idx_[t1], 1), V(idx_[t1], 2));
//        printf("v %f %f %f\n", V(idx_[t2], 0), V(idx_[t2], 1), V(idx_[t2], 2));
//    }
}

// to merge!
void PolylineBranch::project_to_branch_test(const Eigen::Vector3d& p,
                                       int t0, int t1, int t2,
                                       const Eigen::MatrixX3d& V,
                                       Eigen::Vector3d& res_p, int& res_t){
    // printf("%d %d %d\n", t0, t1, t2);
    if (t0 < 0) {  // *starting* vertex of a loop
        t1 = 0;
        t2 = N_-2;
    }
    if (t1 < 0 && t2 >= 0) {  // repair t1
        if (t2 != t0)
            t1 = t2>t0?0:N_-2;
        else { // can do nothing, use the whole branch
            t1 = 0;
            t2 = N_-2;
        }
    }
    if (t1 >= 0 && t2 < 0) {  // repair t2
        if (t1 != t0)
            t2 = t1>t0?0:N_-2;
        else { // can do nothing, use the whole branch
            t1 = 0;
            t2 = N_-2;
        }
    }
    if (t1 < 0 && t2 < 0) {
        t1 = 0;
        t2 = N_-2;
    }

    int lower = t1;
    int upper = t2;
    if (upper < lower) {
        lower = t2;
        upper = t1;
    }
    // printf("%d %d %d (%d -- %d)\n", t0, t1, t2, lower, upper);

    int N = upper - lower + 1; // t0~t0+1, ..., t1~t1+1
    Eigen::MatrixX3d projections(N, 3);
    Eigen::VectorXd dist(N);
    Eigen::Vector3d proj;
    for (int i = 0; i < N; ++i) {
        project_to_line_segment(p, V.row(idx_[lower+i]), V.row(idx_[lower+i+1]), proj);
        projections.row(i) = proj;
        dist(i) = (proj-p).norm();
    }

    double min_dist = std::numeric_limits<double>::max();
    for (int i = 0; i < N; ++i) {
        if (dist(i) < min_dist) {
            min_dist = dist(i);
            res_t = i;
        }
    }
    res_p = projections.row(res_t);
    res_t += lower;
}

void PolylineBranch::project_to_line_segment(const Eigen::Vector3d& P,
                                             const Eigen::Vector3d& A,
                                             const Eigen::Vector3d& B,
                                             Eigen::Vector3d& res) {
    double t = (P-A).dot(B-A)/(B-A).dot(B-A);
    if (t < 0)
        t = 0;
    if (t > 1)
        t = 1;
    res = (1-t) * A + t * B;
}
