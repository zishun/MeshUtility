/*===========================================================================*\
 * Incremental Remeshing
 *
 * The method is described in section 6.5.3 (p.100--p.102) of
 * Mario Botsch, Leif Kobbelt, Mark Pauly, Pierre Alliez, Bruno Lévy. 2010.
 * Polygon Mesh Processing. A K Peters, Ltd.
 *
 * Zishun Liu <liuzishun@gmail.com>
 * July 27, 2021.
\*===========================================================================*/

#ifndef INCREMENTAL_REMESHER_H
#define INCREMENTAL_REMESHER_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>

#include "igl/AABB.h"
#include "PolylineNetwork.h"


class IncrementalRemesher{

    struct EigenTraits : OpenMesh::DefaultTraits {
        using Point = Eigen::Vector3d;
        using Normal = Eigen::Vector3d;
    };
    using EigenTriMesh = OpenMesh::TriMesh_ArrayKernelT<EigenTraits>;

public:

    IncrementalRemesher(){}
    ~IncrementalRemesher() {}

    // mesh IO
    void init_mesh(const Eigen::MatrixXd& _V,
                   const Eigen::MatrixXi& _F);
    void get_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> get_mesh();
    void set_mesh(EigenTriMesh& mesh, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
    void get_mesh(EigenTriMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    void update_points(const Eigen::MatrixXd& V);
    // void find_obtuse_triangle()

    // feature stuff
    void set_features(const Eigen::MatrixXi& E);
    void set_feature_vertices(const Eigen::VectorXi& V);
    void set_features_by_dihedral_angle(const double angle);
    void set_feature_vertices_by_angle(const double angle);
    void collect_features();
    Eigen::MatrixXi get_features();
    Eigen::MatrixXi get_general_features();
    Eigen::VectorXi get_feature_vertices();

    // main steps
    void remesh(const double _targetEdgeLength, const int _iters = 10);
    void split_long_edges(const double _maxEdgeLength );
    void collapse_short_edges(const double _minEdgeLength, const double _maxEdgeLength );
    void equalize_valences();
    void tangential_relaxation();
    void project_to_surface();

private:
    inline int target_valence(const EigenTriMesh::VertexHandle& _vh );
    inline bool is_feature(const EigenTriMesh::EdgeHandle& _eh);
    inline bool is_general_feature(const EigenTriMesh::EdgeHandle& _eh);
    inline bool is_on_feature(const EigenTriMesh::VertexHandle& _vh);
    inline bool is_strong_feature(const EigenTriMesh::VertexHandle& _vh);
    inline bool is_on_general_feature(const EigenTriMesh::VertexHandle& _vh);

private:
    igl::AABB<Eigen::MatrixXd,3>  tree_;
    EigenTriMesh                  mesh_;
    Eigen::MatrixXd               Vin_;
    Eigen::MatrixXi               Fin_;
    bool                          has_feature_edge_ = false;
    bool                          feature_collected_ = true;
    PolylineNetwork               feature_net_;
};

#endif // INCREMENTAL_REMESHER_H defined
