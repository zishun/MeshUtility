// This file is part of MeshUtility, a collection of python mesh processing
// utilities.
//
// Copyright (C) 2021 Zishun Liu <liuzishun@gmail.com>
//
// This Source Code Form is subject to the terms of the 3-Clause BSD License.
// If a copy of the BSD-3-Clause was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/BSD-3-Clause.
#include "mmp_wrapper.h"

mmp_wrapper::mmp_wrapper(){}
mmp_wrapper::~mmp_wrapper(){}


void mmp_wrapper::find_path(const Eigen::MatrixXd& V,
                            const Eigen::MatrixXi& F,
                            int src,
                            int dst,
                            Eigen::MatrixXi& path_edge,
                            Eigen::VectorXd& path_ratio)
{
    std::vector<double> points;
    std::vector<unsigned> faces;
    Eigen::VectorXd D;
    VF_eigenmat_to_stlvector(V, F, points, faces);
    std::vector<geodesic::SurfacePoint> path; // to matrixXd
    find_path(points, faces,
              src, dst,
              path,
              path_edge,
              path_ratio);

    //print_info_about_path(path);
    //    for(unsigned i = 0; i<path.size(); ++i)
    //    {
    //        geodesic::SurfacePoint& s = path[i];
    //        std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
    //    }

    //// output
    // write_path_obj(output_path, path);
    // Eigen::write_binary(output_array, path_info);
}

bool mmp_wrapper::VF_eigenmat_to_stlvector(const Eigen::MatrixXd& V,
                              const Eigen::MatrixXi& F,
                              std::vector<double>& points,
                              std::vector<unsigned>& faces)
{
    points.clear();
    faces.clear();
    points.reserve(V.rows()*3);
    faces.reserve(F.rows()*3);
    for (Eigen::Index i = 0; i < V.rows(); ++i)
    {
        points.push_back(V(i,0));
        points.push_back(V(i,1));
        points.push_back(V(i,2));
    }
    for (Eigen::Index i = 0; i < F.rows(); ++i)
    {
        faces.push_back(F(i,0));
        faces.push_back(F(i,1));
        faces.push_back(F(i,2));
    }

    return true;
}

double mmp_wrapper::point_ratio_on_edge(geodesic::Point3D p0, geodesic::Point3D p1,
                           geodesic::Point3D p)
{
    double d0 = p0.distance(&p);
    double d1 = p1.distance(&p);
    return d0/(d0+d1);
}

void mmp_wrapper::find_path(const std::vector<double>& points,
                            const std::vector<unsigned>& faces,
                            int src, int dst,
                            std::vector<geodesic::SurfacePoint>& path,
                            Eigen::MatrixXi& path_edge,
                            Eigen::VectorXd& path_ratio)
{
    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);

    geodesic::SurfacePoint source(&mesh.vertices()[src]);
    geodesic::SurfacePoint target(&mesh.vertices()[dst]);

    geodesic::GeodesicAlgorithmExact algorithm(&mesh);

    algorithm.geodesic(source, target, path); //find a single source-target path

    path_edge.resize(path.size(), 2);
    path_ratio.resize(path.size());
    for (size_t i = 0; i < path.size(); ++i)
    {
        size_t j = path.size()-i-1; // revert the order
        geodesic::SurfacePoint& p = path[j];
        switch (p.type()) {
        case geodesic::PointType::VERTEX:
            //std::cout << i << ": " << p.base_element()->id() << std::endl;
            path_edge(i,0) = p.base_element()->id();
            path_edge(i,1) = p.base_element()->id();
            path_ratio(i) = 0.0;
            break;
        case geodesic::PointType::EDGE:
        {
            geodesic::base_pointer e = p.base_element();
            geodesic::vertex_pointer v0 = e->adjacent_vertices()[0];
            geodesic::vertex_pointer v1 = e->adjacent_vertices()[1];
            //std::cout << i << ": " << v0->id() << " " << v1->id() << std::endl;
            path_edge(i,0) = v0->id();
            path_edge(i,1) = v1->id();
            path_ratio(i) = point_ratio_on_edge(v0, v1, p);
            break;
        }
        default:
            std::cout << j << ": unknown point type: " << p.type() << std::endl;
            break;
        }
    }
}

double mmp_wrapper::construct_mesh_for_edge_based_geodesic(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXi& source_indices,
        std::vector<double>& points,
        std::vector<unsigned>& faces,
        Eigen::MatrixXi& inserted_vertex_info)
{
    VF_eigenmat_to_stlvector(V, F, points, faces);
    // V, F will not be used from now on.

    // check if 3 source vertices are on the same face.
    Eigen::VectorXi source_vertex_mask;
    source_vertex_mask.setZero(points.size()/3);
    for (Eigen::Index i = 0; i < source_indices.size(); ++i)
    {
        source_vertex_mask[source_indices[i]] = 1;
    }
    for (size_t f = 0; f < faces.size()/3; ++f)
    {
        if (source_vertex_mask[f*3] == 1 && source_vertex_mask[f*3+1] == 1 &&
                source_vertex_mask[f*3+2] == 1)
        {
            // subdivide the face
            //              v0
            //             /|\
            //            / | -\
            //           /  |   -\
            //          /   +     -\
            //         /   /v3\     -\
            //        /  /-    --\    -\
            //       / /-         ---\  -\
            //      //-               --\ -\
            //  v1 /------------------------\  v2
            //
            // generated with https://textik.com/

            int v0 = faces[f*3];
            int v1 = faces[f*3+1];
            int v2 = faces[f*3+2];
            int k = points.size()/3;
            points.push_back((points[v0*3] + points[v1*3] + points[v2*3])/3);
            points.push_back((points[v0*3+1] + points[v1*3+1] + points[v2*3+1])/3);
            points.push_back((points[v0*3+2] + points[v1*3+2] + points[v2*3+2])/3);
            faces[f*3+2] = k;
            faces.push_back(v1); faces.push_back(v2); faces.push_back(k);
            faces.push_back(v2); faces.push_back(v0); faces.push_back(k);
        }
    }

    // construct OpenMesh
    OMesh mesh;
    int nV = points.size()/3;
    std::vector<OMesh::VertexHandle> vhandles(nV);
    for (int i = 0; i < nV; ++i)
    {
        vhandles[i] = mesh.add_vertex(OMesh::Point(points[i*3], points[i*3+1], points[i*3+2]));
    }

    std::vector<OMesh::VertexHandle> face_vhandles(3);
    int nF = faces.size()/3;
    for (int i = 0; i < nF; ++i)
    {
        face_vhandles[0] = vhandles[faces[i*3]];
        face_vhandles[1] = vhandles[faces[i*3+1]];
        face_vhandles[2] = vhandles[faces[i*3+2]];
        mesh.add_face(face_vhandles);
    }

    std::vector<int> inserted_vertex;
    std::vector<int> edge_0;
    std::vector<int> edge_1;

    // find edges
    float max_source_edge_length = 0.0;
    for (int i = 0; i < source_indices.size()-1; ++i)
    {
        OMesh::VertexHandle vh0 = mesh.vertex_handle(source_indices[i]);
        OMesh::VertexHandle vh1 = mesh.vertex_handle(source_indices[i+1]);
        OMesh::HalfedgeHandle he = mesh.find_halfedge(vh0, vh1);
        if (!he.is_valid())
        {
            std::cout << "can not find edge: " << source_indices[i] << "->" << source_indices[i+1] << std::endl;
        }
        else
        {
//         2
//        /|`.
//       / |  `.
//      /  |    `.
//    0 -- k ----- 1
//      \  |    .'
//       \ |  .'
//        \|.'
//         3

            float l = (mesh.point(vh0) - mesh.point(vh1)).norm();
            max_source_edge_length = (l>max_source_edge_length)?l:max_source_edge_length;

            // insert vertex
            points.push_back(0.0);
            points.push_back(0.0);
            points.push_back(0.0);
            int k = points.size()/3-1;
            inserted_vertex.push_back(k);
            edge_0.push_back(vh0.idx());
            edge_1.push_back(vh1.idx());

            // find neighboring faces, and split them
            if (!mesh.is_boundary(he))
            {
                OMesh::FaceHandle fh = mesh.face_handle(he);
                OMesh::HalfedgeHandle nhe = mesh.next_halfedge_handle(he);
                OMesh::VertexHandle vh2 = mesh.to_vertex_handle(nhe);
                int fid = fh.idx();
                faces[3*fid]   = vh1.idx();
                faces[3*fid+1] = vh2.idx();
                faces[3*fid+2] = k;
                faces.push_back(vh2.idx());
                faces.push_back(vh0.idx());
                faces.push_back(k);
            }
            he = mesh.opposite_halfedge_handle(he);
            if (!mesh.is_boundary(he))
            {
                OMesh::FaceHandle fh = mesh.face_handle(he);
                OMesh::HalfedgeHandle nhe = mesh.next_halfedge_handle(he);
                OMesh::VertexHandle vh2 = mesh.to_vertex_handle(nhe);
                int fid = fh.idx();
                faces[3*fid]   = vh0.idx();
                faces[3*fid+1] = vh2.idx();
                faces[3*fid+2] = k;
                faces.push_back(vh2.idx());
                faces.push_back(vh1.idx());
                faces.push_back(k);
            }
        }
    }
    // output: inserted vertex index, two end vertex index
    inserted_vertex_info.resize(inserted_vertex.size(), 3);
    for (size_t i = 0; i < inserted_vertex.size(); ++i)
    {
        inserted_vertex_info(i,0) = inserted_vertex[i];
        inserted_vertex_info(i,1) = edge_0[i];
        inserted_vertex_info(i,2) = edge_1[i];
    }

    // std::cout << "max source edge length: " << max_source_edge_length << std::endl;
    return max_source_edge_length;
}

void mmp_wrapper::edge_sourced_geodesic(std::vector<double>& points,
                                        const std::vector<unsigned>& faces,
                                        int num_dst,
                                        int segments,
                                        Eigen::MatrixXi& inserted_vertex_info,
                                        Eigen::VectorXd& D,
                                        double max_propagation_distance)
{
    geodesic::Mesh mesh;
    D.setConstant(num_dst, std::numeric_limits<double>::max());

    double step = 1.0/segments;
    // We do not compute for #edges * #segments times. Only #segments.
    for (double t = step*0.5; t < 1.0; t += step)
    {
        // update virtual vertex positions
        for (Eigen::Index i = 0; i < inserted_vertex_info.rows(); ++i)
        {
            int vid = inserted_vertex_info(i,0);
            int v1 = inserted_vertex_info(i,1);
            int v2 = inserted_vertex_info(i,2);
            points[vid*3  ] = points[v1*3  ]*t + points[v2*3  ]*(1.0-t);
            points[vid*3+1] = points[v1*3+1]*t + points[v2*3+1]*(1.0-t);
            points[vid*3+2] = points[v1*3+2]*t + points[v2*3+2]*(1.0-t);
        }

        mesh.initialize_mesh_data(points, faces);
        geodesic::GeodesicAlgorithmExact algorithm(&mesh);	//create exact algorithm for the mesh
        std::vector<geodesic::SurfacePoint> all_sources;
        all_sources.reserve(inserted_vertex_info.rows());
        for (Eigen::Index i = 0; i < inserted_vertex_info.rows(); ++i)
        {
            all_sources.push_back(&mesh.vertices()[inserted_vertex_info(i,0)]);
        }

        // still produces 1e100 somewhere!
        algorithm.propagate(all_sources, max_propagation_distance);

        //for(unsigned i=0; i<mesh.vertices().size(); ++i)
        for(unsigned i=0; i<num_dst; ++i)
        {
            geodesic::SurfacePoint p(&mesh.vertices()[i]);

            double distance;
            //unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source
            algorithm.best_source(p,distance);
            //std::cout << distance << " ";		//print geodesic distance for every vertex
            if (distance < D[i])
                D[i] = distance;
        }
        //std::cout << std::endl;
    }
}

void mmp_wrapper::vertex_sourced_geodesic(const std::vector<double>& points,
                                          const std::vector<unsigned>& faces,
                                          const Eigen::VectorXi& source_indices,
                                          Eigen::VectorXd& D,
                                          double max_propagation_distance)
{
    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    geodesic::GeodesicAlgorithmExact algorithm(&mesh);	// create exact algorithm for the mesh
    std::vector<geodesic::SurfacePoint> all_sources;
    all_sources.reserve(source_indices.size());
    for (size_t i = 0; i < source_indices.size(); ++i)
    {
        all_sources.push_back(&mesh.vertices()[source_indices[i]]);
    }

    algorithm.propagate(all_sources, max_propagation_distance);
    D.resize(mesh.vertices().size());
    for(unsigned i=0; i<mesh.vertices().size(); ++i)
    {
        geodesic::SurfacePoint p(&mesh.vertices()[i]);

        double distance;
        //unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source
        algorithm.best_source(p, distance);
        //std::cout << distance << " ";		//print geodesic distance for every vertex
        D[i] = distance;
    }
    //std::cout << std::endl;
}


void mmp_wrapper::distance_field(const Eigen::MatrixXd& V,
                                 const Eigen::MatrixXi& F,
                                 const Eigen::VectorXi& source_indices,
                                 Eigen::VectorXd& D,
                                 double edge_split,
                                 double max_propagation_distance)
{
    if (max_propagation_distance < 0.0)
        max_propagation_distance = std::numeric_limits<double>::max();	// cover the whole mesh

    std::vector<double> points;
    std::vector<unsigned> faces;
    Eigen::MatrixXi inserted_vertex_info;
    if (edge_split > 0.0)
    {
        double max_length =
        construct_mesh_for_edge_based_geodesic(V, F, source_indices,
                                               points, faces, inserted_vertex_info);
        int segments = std::ceil(max_length/edge_split);
        edge_sourced_geodesic(points, faces,
                              V.rows(), // int num_dst,
                              segments,
                              inserted_vertex_info,
                              D,
                              max_propagation_distance);
        // results on sources should be 0.
        for (size_t i = 0; i < source_indices.size(); ++i)
        {
            D[source_indices[i]] = 0.0;
        }
    }
    else
    {
        VF_eigenmat_to_stlvector(V, F, points, faces);
        vertex_sourced_geodesic(points, faces,
                                source_indices,
                                D,
                                max_propagation_distance);
    }
}
