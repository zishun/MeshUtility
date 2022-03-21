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

#ifndef MESHUTILITY_FASTMARCHING_H
#define MESHUTILITY_FASTMARCHING_H

#include <vector>
#include <queue>

#include <Eigen/Dense>
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "OpenMesh/Core/Geometry/EigenVectorT.hh"


class fastmarching 
{
    struct EigenTraits : OpenMesh::DefaultTraits {
        using Point = Eigen::Vector3d;
        //using Normal = Eigen::Vector3d;
        //using TexCoord2D = Eigen::Vector2d;
    };
    using EigenTriMesh = OpenMesh::TriMesh_ArrayKernelT<EigenTraits>;

public:
    fastmarching(const Eigen::MatrixX3d& V, 
                 const Eigen::MatrixX3i& F);
    ~fastmarching(){}

    void set_per_face_gradient_magnitude(Eigen::VectorXd& magnitudes);
    void run(const Eigen::VectorXi& sources, 
             const size_t& n_iter, 
             const double& radio,
             Eigen::VectorXd& dist);

private:
    double update(int & d, const EigenTriMesh::HalfedgeHandle & he, Eigen::Vector3d & vx);
    double planar_update(int & d, Eigen::Matrix3Xd & X, int * x, Eigen::Vector3d & vx, int face_id);

private:
    EigenTriMesh    mesh_;
    double*         per_face_gradient_magnitude_ = nullptr;
    double*         dist_ = nullptr;

};

#endif  // MESHUTILITY_FASTMARCHING_H
