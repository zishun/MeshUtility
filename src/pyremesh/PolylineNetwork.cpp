#include "PolylineNetwork.h"

#include <iostream>
#include <fstream>
#include <string>

void PolylineNetwork::construct_network(const Eigen::MatrixX3d& V,
                                        const Eigen::MatrixXi& E,
                                        const Eigen::VectorXi& S,
                                        Eigen::MatrixXi& Pos) {
    //V_ = V;

    SpMat network(V.rows(), V.rows());
    std::vector<T> triplets;
    triplets.reserve(E.rows()*2);
    for (int i = 0; i < E.rows(); ++i) {
        triplets.push_back(T(E(i,0), E(i,1), 1));
        triplets.push_back(T(E(i,1), E(i,0), 1));
    }
    network.setFromTriplets(triplets.begin(), triplets.end());

    // S to a mask
    std::vector<bool> is_S(V.rows(), false);
    for (int i = 0; i < S.size(); ++i)
        is_S[S[i]] = true;

    // extract branches
    /// branches containing S vertices
    while (true) {
        std::vector<int> idx;

        // find an unvisited edge from a S vertex
        for (int i = 0; i < S.size(); ++i) {
            bool found = false;
            for (SpMat::InnerIterator it(network, S[i]); it; ++it)
            {
                if (it.value() == 1) {
                    idx.push_back(it.row());
                    idx.push_back(it.col());
                    network.coeffRef(it.row(), it.col()) = 2;
                    network.coeffRef(it.col(), it.row()) = 2;
                    found = true;
                    break;
                }
            }
            if (found)
                break;
        }

        if (idx.size() == 0)
            break;

        // extend idx to reach S again
        int end = idx.back();
        while (!is_S[end]) {
            for (SpMat::InnerIterator it(network, end); it; ++it)
            {
                if (it.value() == 1) {
                    idx.push_back(it.col());
                    network.coeffRef(it.row(), it.col()) = 2;
                    network.coeffRef(it.col(), it.row()) = 2;
                    end = it.col();
                    break;
                }
            }
        }

        // add the newly found branch
        branches_.push_back(std::unique_ptr<PolylineBranch>(new PolylineBranch(idx)));
    }

    /// branches containing no S vertices.
    while (true) {
        std::vector<int> idx;

        for (int i = 0; i < network.outerSize(); ++i) {
            bool found = false;
            for (SpMat::InnerIterator it(network, i); it; ++it)
            {
                if (it.value() == 1) {
                    idx.push_back(i);
                    idx.push_back(it.col());
                    network.coeffRef(it.row(), it.col()) = 2;
                    network.coeffRef(it.col(), it.row()) = 2;
                    found = true;
                    break;
                }
            }
            if (found)
                break;
        }

        if (idx.size() == 0)
            break;

        // extend idx
        bool extendable = true;
        int end = idx.back();
        while (extendable) {
            extendable = false;
            for (SpMat::InnerIterator it(network, end); it; ++it)
            {
                if (it.value() == 1) {
                    idx.push_back(it.col());
                    network.coeffRef(it.row(), it.col()) = 2;
                    network.coeffRef(it.col(), it.row()) = 2;
                    end = it.col();
                    extendable = true;
                    break;
                }
            }
        }

        // add the newly found branch
        branches_.push_back(std::unique_ptr<PolylineBranch>(new PolylineBranch(idx)));
    }

    Pos.setConstant(V.rows(), 2, -1);
    for (int i = 0; i < branches_.size(); ++i) {
        // skip two ends. They do not belong to any unique branch.
        for (int j = 1; j < branches_[i]->N_-1; ++j) {
            Pos(branches_[i]->idx_[j], 0) = i;
            Pos(branches_[i]->idx_[j], 1) = j;
        }
    }

    // debug
//    for (int i = 0; i < branches_.size(); ++i) {
//        std::string str = std::to_string(i);
//        std::ofstream out(str+".obj");
//        // skip two ends. They do not belong to any unique branch.
//        for (int j = 0; j < branches_[i]->N_; ++j) {
//            out << "v " <<
//                   V(branches_[i]->idx_[j], 0) << " " <<
//                   V(branches_[i]->idx_[j], 1) << " " <<
//                   V(branches_[i]->idx_[j], 2) << std::endl;
//        }
//        for (int j = 0; j < branches_[i]->N_-1; ++j) {
//            out << "l " << j+1 << " " << j+2  << std::endl;
//        }
//        out.close();
//    }
}

void PolylineNetwork::project_to_network(const Eigen::Vector3d& p, const int branch,
                                         int t0, int t1, int t2,
                                         const Eigen::MatrixX3d& V,
                                         Eigen::Vector3d& res_p, int& res_t) {
    branches_[branch]->project_to_branch(p, t0, t1, t2, V, res_p, res_t);
}

int PolylineNetwork::find_branch(int v0, int v1) {
    // TODO: what if there are multiple matchings?
    for (int i = 0; i < branches_.size(); ++i) {
        if (branches_[i]->idx_[0] == v0 && branches_[i]->idx_.back() == v1)
            return i;
        if (branches_[i]->idx_[0] == v1 && branches_[i]->idx_.back() == v0)
            return i;
    }
    return -1;  // TODO: are you going to handle this explicitly?
}

int PolylineNetwork::find_segment(int branch, int v) {
    // TODO: what if this is a loop and v is the *starting* vertex?
    if (branches_[branch]->idx_[0] == v)
        return 0;
    else
        return branches_[branch]->N_ - 2;
}

void PolylineNetwork::print_network() {
    for (size_t i = 0; i < branches_.size(); ++i) {
        for (size_t j = 0; j < branches_[i]->idx_.size(); ++j) {
            std::cout << branches_[i]->idx_[j] << " ";
        }
        std::cout << std::endl;
    }
}
