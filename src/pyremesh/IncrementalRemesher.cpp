#include "IncrementalRemesher.h"

void IncrementalRemesher::init_mesh(const Eigen::MatrixXd& _V,
                                    const Eigen::MatrixXi& _F)
{
    tree_.init(_V, _F);
    Vin_ = _V;
    Fin_ = _F;
    set_mesh(mesh_, _V, _F);

    mesh_.request_edge_status();
    mesh_.request_vertex_status();
    mesh_.request_face_status();
    mesh_.request_face_normals();
    mesh_.request_vertex_normals();
    mesh_.update_normals();
}

void IncrementalRemesher::set_mesh(EigenTriMesh& _mesh, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    _mesh.clean();
    std::vector<EigenTriMesh::VertexHandle> vhandles(V.rows());
    for (Eigen::Index i = 0; i < V.rows(); ++i)
    {
        vhandles[i] = _mesh.add_vertex(V.row(i));
    }

    std::vector<EigenTriMesh::VertexHandle> face_vhandles(3);
    for (Eigen::Index i = 0; i < F.rows(); ++i)
    {
        face_vhandles[0] = vhandles[F(i,0)];
        face_vhandles[1] = vhandles[F(i,1)];
        face_vhandles[2] = vhandles[F(i,2)];
        _mesh.add_face(face_vhandles);
    }

    for (auto vh : _mesh.all_vertices()) {
        if (_mesh.is_boundary(vh)) {
            has_feature_edge_ = true;
            feature_collected_ = false;
            break;
        }
    }
}

void IncrementalRemesher::get_mesh(EigenTriMesh& _mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    V.resize(_mesh.n_vertices(), 3);
    F.resize(_mesh.n_faces(), 3);
    for (auto vh : _mesh.all_vertices()) {
        EigenTriMesh::Point p = _mesh.point(vh);
        V(vh.idx(), 0) = p[0];
        V(vh.idx(), 1) = p[1];
        V(vh.idx(), 2) = p[2];
    }
    for (auto fh : _mesh.all_faces()) {
        auto fv_it = _mesh.fv_iter(fh);
        F(fh.idx(), 0) = fv_it->idx(); ++fv_it;
        F(fh.idx(), 1) = fv_it->idx(); ++fv_it;
        F(fh.idx(), 2) = fv_it->idx();
    }
}

void IncrementalRemesher::get_mesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    get_mesh(mesh_, V, F);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> IncrementalRemesher::get_mesh()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    get_mesh(mesh_, V, F);
    auto res = std::make_tuple(V, F);
    return res;
}

void IncrementalRemesher::update_points(const Eigen::MatrixXd& V)
{
    if (mesh_.n_vertices() != V.rows()) {
        std::cerr << "Cannot update_points, incorrect number of vertices!" << std::endl;
        return;
    }

    for (auto vh : mesh_.all_vertices()) {
        mesh_.set_point(vh, V.row(vh.idx()));
    }
}

void IncrementalRemesher::set_features(const Eigen::MatrixXi& _E)
{
    feature_collected_ = false;

    for (int i = 0; i < _E.rows(); ++i)
    {
        EigenTriMesh::VertexHandle vh0 = mesh_.vertex_handle(_E(i,0));
        EigenTriMesh::VertexHandle vh1 = mesh_.vertex_handle(_E(i,1));
        EigenTriMesh::HalfedgeHandle heh = mesh_.find_halfedge(vh0, vh1);
        if (heh.is_valid())
        {
            EigenTriMesh::EdgeHandle eh = mesh_.edge_handle(heh);
            mesh_.status(eh).set_feature(true);
        }
        else
        {
            std::cout << "Invalid feature edge: " << _E.row(i) << std::endl;
        }
    }

    if (_E.rows() > 0)
        has_feature_edge_ = true;

    // find strong feature vertices automatically by:
    set_feature_vertices_by_angle(180.0);
}

void IncrementalRemesher::set_feature_vertices(const Eigen::VectorXi& V)
{
    feature_collected_ = false;

    for (int i = 0; i < V.size(); ++i)
    {
        EigenTriMesh::VertexHandle vh = mesh_.vertex_handle(V[i]);
        mesh_.status(vh).set_feature(true);
    }
}

void IncrementalRemesher::set_features_by_dihedral_angle(const double angle)
{
    feature_collected_ = false;

    const double cos_tol = std::cos(angle/180.0*M_PI);

    for (auto eh : mesh_.all_edges()) {
        if (mesh_.is_boundary(eh))
            continue;
        const EigenTriMesh::HalfedgeHandle & hh0 = mesh_.halfedge_handle(eh, 0);
        const EigenTriMesh::HalfedgeHandle & hh1 = mesh_.halfedge_handle(eh, 1);
        const EigenTriMesh::FaceHandle & fh0 = mesh_.face_handle(hh0);
        const EigenTriMesh::FaceHandle & fh1 = mesh_.face_handle(hh1);
        Eigen::Vector3d n0 = mesh_.normal(fh0);
        Eigen::Vector3d n1 = mesh_.normal(fh1);
        if (n0.dot(n1) < cos_tol) {
            mesh_.status( eh ).set_feature( true );
            has_feature_edge_ = true;
        }
    }
}

void IncrementalRemesher::set_feature_vertices_by_angle(const double angle)
{
    feature_collected_ = false;

    const double cos_tol = std::cos(angle/180.0*M_PI);

    for (auto vh : mesh_.all_vertices()) {
        int cnt = 0;
        Eigen::Vector3d vec0;
        Eigen::Vector3d vec1;
        EigenTriMesh::VOHIter vh_it;
        for (vh_it = mesh_.voh_iter(vh); vh_it.is_valid(); ++vh_it) {
            const EigenTriMesh::EdgeHandle & eh = mesh_.edge_handle( *vh_it );
            if (is_general_feature(eh)) {
                const EigenTriMesh::VertexHandle & v1 = mesh_.to_vertex_handle(*vh_it);
                EigenTriMesh::Point vec = mesh_.point(v1) - mesh_.point(vh);

                if (cnt == 0) {
                    vec0 = vec / vec.norm();
                } else {
                    vec1 = vec / vec.norm();
                }
                ++cnt;
            }
        }
        if (cnt > 0) {
            if (cnt == 2) {
                if (vec0.dot(vec1) > -cos_tol) {
                    mesh_.status( vh ).set_feature( true );
                }
            }else{ // i.e. corner vertices on page.103.
                mesh_.status( vh ).set_feature( true );
            }
        }
    }
}

void IncrementalRemesher::collect_features() {
    if (!has_feature_edge_)
        return;

    // make sure joints are marked as strong features.
    // find them automatically by:
    set_feature_vertices_by_angle(180.0);

    Eigen::MatrixXi E = get_general_features();
    Eigen::VectorXi S = get_feature_vertices();
    Eigen::MatrixXi Pos;
    feature_net_.construct_network(Vin_, E, S, Pos);

    // set Pos;
    OpenMesh::VPropHandleT< Eigen::Vector2i > feature_pos;
    if ( !mesh_.get_property_handle(feature_pos, "Feature Position") )
        mesh_.add_property(feature_pos, "Feature Position" );
    for (auto vh : mesh_.all_vertices()) {
        mesh_.property(feature_pos, vh) = Pos.row(vh.idx());
    }

//    for (auto vh : mesh_.all_vertices()) {
//        if (is_on_general_feature(vh))
//            std::cout << mesh_.property(feature_pos, vh) << std::endl;
//    }
    feature_collected_ = true;
}

Eigen::MatrixXi IncrementalRemesher::get_features()
{
    std::vector<int> idx;
    for (auto eh : mesh_.all_edges()) {
        if (is_feature(eh)) {
            const EigenTriMesh::HalfedgeHandle & hh = mesh_.halfedge_handle(eh, 0);
            const EigenTriMesh::VertexHandle & v0 = mesh_.from_vertex_handle(hh);
            const EigenTriMesh::VertexHandle & v1 = mesh_.to_vertex_handle(hh);
            idx.push_back(v0.idx());
            idx.push_back(v1.idx());
        }
    }
    Eigen::MatrixXi E = Eigen::Map<Eigen::MatrixXi>(idx.data(), 2, idx.size()/2);
    return E.transpose();
}

Eigen::MatrixXi IncrementalRemesher::get_general_features()
{
    std::vector<int> idx;
    for (auto eh : mesh_.all_edges()) {
        if (is_general_feature(eh)) {
            const EigenTriMesh::HalfedgeHandle & hh = mesh_.halfedge_handle(eh, 0);
            const EigenTriMesh::VertexHandle & v0 = mesh_.from_vertex_handle(hh);
            const EigenTriMesh::VertexHandle & v1 = mesh_.to_vertex_handle(hh);
            idx.push_back(v0.idx());
            idx.push_back(v1.idx());
        }
    }
    Eigen::MatrixXi E = Eigen::Map<Eigen::MatrixXi>(idx.data(), 2, idx.size()/2);
    return E.transpose();
}

Eigen::VectorXi IncrementalRemesher::get_feature_vertices()
{
    std::vector<int> idx;
    for (auto vh : mesh_.all_vertices()) {
        if (is_strong_feature(vh))
            idx.push_back(vh.idx());
    }

    Eigen::VectorXi V = Eigen::Map<Eigen::VectorXi>(idx.data(), idx.size(), 1);
    return V;
}

/// main loop. pseudocode at p.100 middle.
void IncrementalRemesher::remesh(const double _targetEdgeLength,
                                 const int _iters)
{
    if (!feature_collected_)
        collect_features();

    const double low  = (4.0 / 5.0) * _targetEdgeLength;
    const double high = (4.0 / 3.0) * _targetEdgeLength;

    for (int i=0; i < _iters; i++){
        split_long_edges(high);
        collapse_short_edges(low, high);
        equalize_valences();
        tangential_relaxation();
        project_to_surface();
    }
}

/// performs edge splits until all edges are shorter than the threshold
/// pseudocode at p.100 bottom.
void IncrementalRemesher::split_long_edges( const double _maxEdgeLength )
{
    if (!feature_collected_)
        collect_features();

    OpenMesh::VPropHandleT< Eigen::Vector2i > feat_pos;
    if (has_feature_edge_)
        mesh_.get_property_handle(feat_pos, "Feature Position");

    EigenTriMesh::EdgeIter e_it;
    EigenTriMesh::EdgeIter e_end = mesh_.edges_end();

    for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it){
        const EigenTriMesh::HalfedgeHandle & hh = mesh_.halfedge_handle( *e_it, 0 );

        const EigenTriMesh::VertexHandle & v0 = mesh_.from_vertex_handle(hh);
        const EigenTriMesh::VertexHandle & v1 = mesh_.to_vertex_handle(hh);

        EigenTriMesh::Point vec = mesh_.point(v1) - mesh_.point(v0);

        if ( vec.norm() > _maxEdgeLength ){
            const EigenTriMesh::Point midPoint = mesh_.point(v0) + ( 0.5 * vec );

            // split at midpoint
            EigenTriMesh::VertexHandle vh = mesh_.add_vertex( midPoint );

            bool had_general_feature = is_general_feature(*e_it); // get it before splitting
            bool had_feature = is_feature(*e_it); // get it before splitting

            mesh_.split(*e_it, vh);

            if ( had_general_feature ){
                EigenTriMesh::VOHIter vh_it;
                if (had_feature) {
                    for (vh_it = mesh_.voh_iter(vh); vh_it.is_valid(); ++vh_it)
                        if ( mesh_.to_vertex_handle(*vh_it) == v0 || mesh_.to_vertex_handle(*vh_it) == v1 )
                            mesh_.status( mesh_.edge_handle( *vh_it ) ).set_feature( true );
                }

                Eigen::Vector2i& pos0 = mesh_.property(feat_pos, v0);
                Eigen::Vector2i& pos1 = mesh_.property(feat_pos, v1);
                int branch = pos0[0];
                int t = pos0[1];
                // TODO: is vertex idx reliable?
                if (branch == -1) {
                    if (pos1[0] == -1) {  // both ends.
                        branch = feature_net_.find_branch(v0.idx(), v1.idx()); // TODO
                        t = 0;
                    } else {
                        branch = pos1[0];  // v0 is an end.
                        t = feature_net_.find_segment(branch, v0.idx()); // TODO
                    }
                } else {
                    if (pos1[0] == -1) {  // v1 is an end.
                        t = feature_net_.find_segment(branch, v1.idx());
                    } else {  // neither is end.
                        t = std::min(pos0[1], pos1[1]);
                    }
                }
                mesh_.property(feat_pos, vh) = Eigen::Vector2i(branch, t);
                //std::cout << "split_long_edges " << branch << " " << t << std::endl;
            }
        }
    }
}

/// collapse edges shorter than minEdgeLength if collapsing doesn't result in new edge longer than maxEdgeLength
/// pseudocode on p.101.
void IncrementalRemesher::collapse_short_edges( const double _minEdgeLength, const double _maxEdgeLength )
{
    if (!feature_collected_)
        collect_features();

    // add checked property, to determine when to exit the loop
    OpenMesh::EPropHandleT< bool > checked;
    if ( !mesh_.get_property_handle(checked, "Checked Property") )
        mesh_.add_property(checked, "Checked Property" );
    for (auto eh : mesh_.all_edges()) {
        mesh_.property(checked, eh) = false;
    }

    // TODO temp
    OpenMesh::VPropHandleT< Eigen::Vector2i > feature_pos;
    mesh_.get_property_handle(feature_pos, "Feature Position");

    bool finished = false;
    while( !finished ){
        finished = true;

        for (auto eh : mesh_.all_edges()) {
            if ( mesh_.property(checked, eh) )
                continue;
            mesh_.property(checked, eh) = true;

            EigenTriMesh::HalfedgeHandle hh = mesh_.halfedge_handle(eh, 0);
            EigenTriMesh::VertexHandle v0 = mesh_.from_vertex_handle(hh);
            EigenTriMesh::VertexHandle v1 = mesh_.to_vertex_handle(hh);

            // check if v0/v1 are feature vertices.
            // TODO: convered all cases?
            if (is_strong_feature(v0) && is_strong_feature(v1))
                continue;
            if (is_strong_feature(v0)) {  // but v1 is not. swap v0/v1.
                hh = mesh_.halfedge_handle(eh, 1);
                v0 = mesh_.from_vertex_handle(hh);
                v1 = mesh_.to_vertex_handle(hh);
            }
            // now v0 is not a strong feature.
            if (is_on_general_feature(v0)) {
                // then only one case is possible: v0v1 is a general feature edge
                // (both feature and on the same feature!)
                if (!(is_on_general_feature(v1) && is_general_feature(eh))) {
                    continue;
                }
            }

            const EigenTriMesh::Point vec = mesh_.point(v1) - mesh_.point(v0);
            const double edgeLength = vec.norm();

            // edge too short but don't try to collapse edges that have length 0
            // skip those already collapsed?
            if ( (edgeLength < _minEdgeLength) && (edgeLength > DBL_EPSILON) ){
            // if ( (edgeLength < _minEdgeLength) ){

                // check if the collapse is ok
                bool collapse_ok = true;

                const EigenTriMesh::Point & B = mesh_.point(v1);

                // allow sliding along feature/boundary edges.
                std::vector<EigenTriMesh::VertexHandle> neighbor_feat;
                if (collapse_ok) {
                    for(EigenTriMesh::VOHIter vh_it = mesh_.voh_iter(v0); vh_it.is_valid(); ++vh_it) {
                        if ( (( B - mesh_.point( mesh_.to_vertex_handle(*vh_it) ) ).norm() > _maxEdgeLength )
                             ){
                            collapse_ok = false;
                            break;
                        }

                        if (is_feature(mesh_.edge_handle( *vh_it ))) {
                            neighbor_feat.push_back(mesh_.to_vertex_handle(*vh_it));
                        }
                    }
                }

//                bool v1_is_feat = is_on_general_feature(v1);
//                if (v1_is_feat)
//                    std::cout << mesh_.property(feature_pos, v1) << std::endl;

                if( collapse_ok && mesh_.is_collapse_ok(hh) ) {
                    mesh_.collapse( hh );

                    for (auto v2 : neighbor_feat) {
                        EigenTriMesh::HalfedgeHandle heh = mesh_.find_halfedge(v1, v2);
                        if (heh.is_valid()) {
                            EigenTriMesh::EdgeHandle eh = mesh_.edge_handle(heh);
                            mesh_.status(eh).set_feature(true);
                        }
                    }

                    finished = false;

//                    if (v1_is_feat)
//                        std::cout << mesh_.property(feature_pos, v1) << std::endl;
                }
            }
        }
    }

    mesh_.remove_property(checked);

    mesh_.garbage_collection();
}

/// pseudocode at p.102 top.
void IncrementalRemesher::equalize_valences()
{
    if (!feature_collected_)
        collect_features();

    EigenTriMesh::EdgeIter e_it;
    EigenTriMesh::EdgeIter e_end = mesh_.edges_end();

    for (e_it = mesh_.edges_sbegin(); e_it != e_end; ++e_it) {
        if (!mesh_.is_flip_ok(*e_it)) continue;
        if (is_feature(*e_it)) continue;

        const EigenTriMesh::HalfedgeHandle & h0 = mesh_.halfedge_handle( *e_it, 0 );
        const EigenTriMesh::HalfedgeHandle & h1 = mesh_.halfedge_handle( *e_it, 1 );

        if (h0.is_valid() && h1.is_valid()) {
            if (mesh_.face_handle(h0).is_valid() && mesh_.face_handle(h1).is_valid()) {
                // vertices of the two triangles adjacent to e_it
                const EigenTriMesh::VertexHandle & a = mesh_.to_vertex_handle(h0);
                const EigenTriMesh::VertexHandle & b = mesh_.to_vertex_handle(h1);
                const EigenTriMesh::VertexHandle & c = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h0));
                const EigenTriMesh::VertexHandle & d = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h1));

                const int deviation_pre =  abs((int)(mesh_.valence(a) - target_valence(a)))
                        +abs((int)(mesh_.valence(b) - target_valence(b)))
                        +abs((int)(mesh_.valence(c) - target_valence(c)))
                        +abs((int)(mesh_.valence(d) - target_valence(d)));
                mesh_.flip(*e_it);

                const int deviation_post = abs((int)(mesh_.valence(a) - target_valence(a)))
                        +abs((int)(mesh_.valence(b) - target_valence(b)))
                        +abs((int)(mesh_.valence(c) - target_valence(c)))
                        +abs((int)(mesh_.valence(d) - target_valence(d)));

                if (deviation_pre <= deviation_post)  // reverse the flipping
                    mesh_.flip(*e_it);
            }
        }
    }
}

/// pseudocode at p.102 bottom.
void IncrementalRemesher::tangential_relaxation()
{
    if (!feature_collected_)
        collect_features();

    mesh_.update_normals();

    // add checked property
    OpenMesh::VPropHandleT< EigenTriMesh::Point > q;
    if ( !mesh_.get_property_handle(q, "q Property") )
        mesh_.add_property(q, "q Property" );

    EigenTriMesh::ConstVertexIter v_it;
    EigenTriMesh::ConstVertexIter v_end = mesh_.vertices_end();

    // smoothing with uniform Laplacian weights (barycenters)
    for (v_it = mesh_.vertices_sbegin(); v_it != v_end; ++v_it){

        if (is_strong_feature(*v_it) || is_on_general_feature(*v_it))
            continue;

        EigenTriMesh::Point tmp(0.0, 0.0, 0.0);
        unsigned int N = 0;

        EigenTriMesh::VVIter vv_it;
        for (vv_it = mesh_.vv_iter(*v_it); vv_it.is_valid(); ++vv_it){
            tmp += mesh_.point(*vv_it);
            N++;
        }

        if (N > 0)
            tmp /= (double) N;

        mesh_.property(q, *v_it) = tmp;
    }

    // move to new position
    OpenMesh::VPropHandleT< Eigen::Vector2i > feature_pos;
    mesh_.get_property_handle(feature_pos, "Feature Position");
    for (v_it = mesh_.vertices_sbegin(); v_it != v_end; ++v_it){
        if (is_strong_feature(*v_it))
            continue;
        if (!is_on_general_feature(*v_it))
            mesh_.set_point(*v_it, mesh_.property(q, *v_it) + (mesh_.normal(*v_it).dot((mesh_.point(*v_it) - mesh_.property(q, *v_it) )) ) * mesh_.normal(*v_it));
        else { // on_general_feature but not strong feature
            EigenTriMesh::VertexHandle neighbors[2];
            int cnt = 0;
            EigenTriMesh::VOHIter vh_it;
            for (vh_it = mesh_.voh_iter(*v_it); vh_it.is_valid(); ++vh_it) {
                const EigenTriMesh::EdgeHandle & eh = mesh_.edge_handle( *vh_it );
                if (is_general_feature(eh)) {
                    EigenTriMesh::VertexHandle v1 = mesh_.to_vertex_handle(*vh_it);
                    neighbors[cnt] = v1;
                    ++cnt;
                }
            }

            // project midpoint to the original feature curve.
            // EigenTriMesh::Point p0 = mesh_.point(*v_it);
            EigenTriMesh::Point p1 = mesh_.point(neighbors[0]);
            EigenTriMesh::Point p2 = mesh_.point(neighbors[1]);
            EigenTriMesh::Point res;
            Eigen::Vector2i& pos0 = mesh_.property(feature_pos, *v_it);
            Eigen::Vector2i& pos1 = mesh_.property(feature_pos, neighbors[0]);
            Eigen::Vector2i& pos2 = mesh_.property(feature_pos, neighbors[1]);
            int branch = pos0[0];
            // usualy v_it have a unique branch id.
            // BUT it may be the *starting* vertex of a loop.
            if (branch == -1)
                branch = pos1[0];
            // printf("branch %d %d %d\n", pos0[0], pos1[0], pos2[0]);
            int t0 = pos0[1];
            int t1 = pos1[1];
            int t2 = pos2[1];
            int res_t;
            feature_net_.project_to_network(0.5*(p1+p2), branch, t0, t1, t2,
                                            Vin_, res, res_t);
            mesh_.set_point(*v_it, res);
            mesh_.property(feature_pos, *v_it) = Eigen::Vector2i(branch, res_t);
        }
    }

    mesh_.remove_property(q);
}

void IncrementalRemesher::project_to_surface()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    get_mesh(mesh_, V, F);

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    tree_.squared_distance(Vin_, Fin_, V, sqrD, I, C);

    for (auto vh : mesh_.all_vertices()) {
        mesh_.set_point(vh, C.row(vh.idx()));
    }
}

/// returns 4 for boundary vertices and 6 otherwise
inline
int IncrementalRemesher::target_valence( const EigenTriMesh::VertexHandle& _vh ){
    return mesh_.is_boundary(_vh)?4:6;
}

inline
bool IncrementalRemesher::is_feature( const EigenTriMesh::EdgeHandle& _eh ){
    return mesh_.status(_eh).feature();
}

inline
bool IncrementalRemesher::is_general_feature( const EigenTriMesh::EdgeHandle& _eh ){
    return is_feature(_eh) || mesh_.is_boundary(_eh);
}

inline
bool IncrementalRemesher::is_on_feature( const EigenTriMesh::VertexHandle& _vh ){

    EigenTriMesh::ConstVertexOHalfedgeIter voh_it;

    for (voh_it = mesh_.voh_iter( _vh ); voh_it.is_valid(); ++voh_it )
        if ( mesh_.status( mesh_.edge_handle( *voh_it ) ).feature() )
            return true;

    return false;
}

inline
bool IncrementalRemesher::is_strong_feature( const EigenTriMesh::VertexHandle& _vh ){
    return mesh_.status(_vh).feature();
}

inline
bool IncrementalRemesher::is_on_general_feature( const EigenTriMesh::VertexHandle& _vh ){
    return is_on_feature(_vh) || mesh_.is_boundary(_vh);
}
