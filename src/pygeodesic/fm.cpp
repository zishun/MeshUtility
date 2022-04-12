// This file is a partial conversion of
// https://github.com/larc/gproshan/blob/main/src/geodesics.cpp
// gproshan: a geometry processing and shape analysis framework
// is released under MIT license.
/*
MIT License

Copyright (c) 2018 Luciano Arnaldo Romero Calla, Lizeth Joseline Fuentes Perez

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "fm.h"


fastmarching::fastmarching(const Eigen::MatrixX3d& V, 
                           const Eigen::MatrixX3i& F)
{
    mesh_.clean();
    std::vector<EigenTriMesh::VertexHandle> vhandles(V.rows());
    for (Eigen::Index i = 0; i < V.rows(); ++i)
    {
        vhandles[i] = mesh_.add_vertex(V.row(i));
    }

    std::vector<EigenTriMesh::VertexHandle> face_vhandles(3);
    for (Eigen::Index i = 0; i < F.rows(); ++i)
    {
        face_vhandles[0] = vhandles[F(i,0)];
        face_vhandles[1] = vhandles[F(i,1)];
        face_vhandles[2] = vhandles[F(i,2)];
        mesh_.add_face(face_vhandles);
    }
}

void fastmarching::set_per_face_gradient_magnitude(Eigen::VectorXd& magnitudes)
{
    if (magnitudes.size() == mesh_.n_faces())
        per_face_gradient_magnitude_ = magnitudes.data();
}

void fastmarching::run(const Eigen::VectorXi& sources, 
                       const size_t& n_iter, 
                       const double& radio,
                       Eigen::VectorXd& dist)
{
    size_t n_vertices = mesh_.n_vertices();
    
    dist.resize(n_vertices);
    dist_ = &dist(0);
    int BLACK = 0, GREEN = 1, RED = 2;
    int * color = new int[n_vertices];

    //#pragma omp parallel for
    for(int v = 0; v < n_vertices; v++)
    {
        color[v] = GREEN;
        dist[v] = std::numeric_limits<double>::max();
    }

    size_t green_count = n_iter ? n_iter : n_vertices;

    std::priority_queue<std::pair<double, size_t>,
            std::vector<std::pair<double, size_t> >,
            std::greater<std::pair<double, size_t> > > Q;

    double dv, dp;
    int dir; // dir propagation
    Eigen::Vector3d vx;

    size_t black_i, v;
    EigenTriMesh::VertexHandle vh;

    //int c = 0;
    //n_sorted = 0;
    for(Eigen::Index i = 0; i < sources.size(); ++i)
    {
        int s = sources[i];
        dist[s] = 0;
        //if(clusters) clusters[s] = ++c;
        color[s] = RED;
        Q.push(std::make_pair(dist[s], s));
    }

    while(green_count-- && !Q.empty())
    {
        while(!Q.empty() && color[Q.top().second] == BLACK)
            Q.pop();

        if(Q.empty()) break;

        black_i = Q.top().second;
        color[black_i] = BLACK;
        Q.pop();
        
        if(dist[black_i] > radio) break;

        //sorted_index[n_sorted++] = black_i;

        // std::vector<size_t> black_link;
        // mesh_.link(black_link, black_i); // TODO
        // for(const int & he: black_link)
        // {
        //     v = mesh->vt(he);
        EigenTriMesh::VertexHandle black_i_vh = mesh_.vertex_handle(black_i);
        for (EigenTriMesh::VertexVertexIter vv_it = mesh_.vv_iter(black_i_vh); vv_it.is_valid(); ++vv_it)
        {
            vh = *vv_it;
            v = vh.idx();

            if(color[v] == GREEN)
                color[v] = RED;

            if(color[v] == RED)
            {
                dv = dist[v];
                //for_star(v_he, mesh, v)
                EigenTriMesh::VertexOHalfedgeIter voh_it = mesh_.voh_iter(vh);
                for(; voh_it.is_valid(); ++voh_it)
                {
                    EigenTriMesh::HalfedgeHandle heh = *voh_it;
                    if (mesh_.is_boundary(heh))
                        continue;
                    dp = update(dir, heh, vx);
                    //dp = update_step(mesh, dist, v_he);
                    if(dp < dv)
                    {
                        dv = dp;
                        
                        //if(clusters)
                        //    clusters[v] = dist[mesh_.vt(prev(v_he))] < dist[mesh_.vt(next(he))] ? clusters[mesh_.vt(prev(he))] : clusters[mesh_.vt(next(he))];
                    }
                }

                if(dv < dist[v])
                    Q.push(std::make_pair(dist[v] = dv, v));
            }
        }
    }

    delete [] color;
}

//d = {NIL, 0, 1} cross edge, next, prev
double fastmarching::update(int & d, const EigenTriMesh::HalfedgeHandle & he, Eigen::Vector3d & vx)
{
    //d = NIL;

    Eigen::Matrix3Xd X(3,2);
    int x[3];

    // x[0] = mesh_.vt(next(he));
    // x[1] = mesh_.vt(prev(he));
    // x[2] = mesh_.vt(he);                //update x[2]
    x[0] = mesh_.to_vertex_handle(he).idx();
    EigenTriMesh::HalfedgeHandle nhe = mesh_.next_halfedge_handle(he);
    x[1] = mesh_.to_vertex_handle(nhe).idx();
    x[2] = mesh_.from_vertex_handle(he).idx();

    vx = mesh_.point(mesh_.vertex_handle(x[2]));
    X.col(0) = mesh_.point(mesh_.vertex_handle(x[0])) - vx;
    X.col(1) = mesh_.point(mesh_.vertex_handle(x[1])) - vx;

    int face = mesh_.face_handle(he).idx();
    return planar_update(d, X, x, vx, face);
}

// Sec 3.1 Planar wavefront approximation (p.5)
// Parallel algorithms for approximation of distance maps on parametric surfaces
// TOG 2008.
// https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.435.6617&rep=rep1&type=pdf
double fastmarching::planar_update(int & d, Eigen::Matrix3Xd & X, int * x, Eigen::Vector3d & vx, int face_id)
{
    Eigen::MatrixXd ones(2,1);
    ones.setOnes(2,1);

    Eigen::MatrixXd Q;
    //if(!inv_sympd(Q, X.transpose() * X))  // Q = (X^T X)^(-1)
    //    return std::numeric_limits<double>::max();
    Q = (X.transpose() * X).inverse();  // Q = (X^T X)^(-1)
    // TODO: check the error.

    Eigen::MatrixXd t(2,1);

    if (per_face_gradient_magnitude_ == nullptr) {
        t(0,0) = dist_[x[0]];
        t(1,0) = dist_[x[1]];
    } else {
        //std::cout << per_face_gradient_magnitude_[face_id] << " ";
        t(0,0) = dist_[x[0]] / per_face_gradient_magnitude_[face_id];
        t(1,0) = dist_[x[1]] / per_face_gradient_magnitude_[face_id];
    }

    // solve eq(10).
    double p;
    Eigen::MatrixXd delta = ones.transpose() * Q * t;
    double dis = (delta * delta - (ones.transpose() * Q * ones) * ((t.transpose() * Q * t)(0) - 1))(0);

    if(dis >= 0)
    {
        p = delta(0) + std::sqrt(dis);
        p /= (ones.transpose() * Q * ones)(0);
    }
    else p = std::numeric_limits<double>::max();
    // end of solving eq(10).

    Eigen::MatrixXd n = X * Q * (t - p * ones);  // eq(8)
    Eigen::MatrixXd cond = Q * X.transpose() * n;  // eq(15)

    Eigen::Vector3d v;

    if(t(0) == std::numeric_limits<double>::max() || t(1) == std::numeric_limits<double>::max() || dis < 0 || (cond(0) >= 0 || cond(1) >= 0))
    {
        // eq(16)
        double dp[2];
        if (per_face_gradient_magnitude_ == nullptr) {
            dp[0] = dist_[x[0]] + X.col(0).norm();
            dp[1] = dist_[x[1]] + X.col(1).norm();
        } else {
            dp[0] = dist_[x[0]] / per_face_gradient_magnitude_[face_id] + X.col(0).norm();
            dp[1] = dist_[x[1]] / per_face_gradient_magnitude_[face_id] + X.col(1).norm();
        }

        d = dp[1] < dp[0];
        v = X.col(d);
        p = dp[d];
    }
    else
    {
        Eigen::MatrixXd A(3,2);
        A.col(0) = -n;
        A.col(1) = X.col(1) - X.col(0);
        Eigen::VectorXd b = -X.col(0);
        Eigen::MatrixXd l = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        v = l(1) * A.col(1) + X.col(0);
    }

    //vx += *((vertex *) v.memptr()); // TODO??

    if (per_face_gradient_magnitude_ != nullptr) {
         p *= per_face_gradient_magnitude_[face_id];
    }

    return p;
}

