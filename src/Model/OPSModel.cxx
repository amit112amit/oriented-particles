#include "OPSModel.h"

namespace OPS{

//! Constructor for BrownOPS
OPSModel::OPSModel(size_t n, double_t &f, RefVXd x, RefVXd g, RefVXd pX):
        _positions(x.data(),3,n), _rotationVectors(&(x(3*n)),3,n),
        _posGradient(g.data(),3,n), _rotGradient(&(g(3*n)),3,n),
        _prevX(pX.data(),3,n), Model(6*n,f,g){

    // Set number of particles
    _N = n;

    // Initialize internal arrays
    _diffNormalRV = std::vector< Matrix3d,
                    Eigen::aligned_allocator<Matrix3d> >(_N, Matrix3d::Zero());
    _first_ring.resize(_N);


    //Convert rotation vectors to point normals
    _normals = Matrix3Xd::Zero(3,_N);
    computeNormals();

    // Construct triangles and edges
    updateTriangles();

    //Set the initial nearest neighbor map
    void *coords = (void*) _positions.data();
    vtkNew<vtkDoubleArray> pointCoords;
    pointCoords->SetVoidArray( coords, 3*_N, 1);
    pointCoords->SetNumberOfComponents(3);

    vtkNew<vtkPoints> points;
    points->SetData( pointCoords );

    // Construct vtkPolyData
    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints( points );

    //Construct the kd-tree
    vtkNew<vtkOctreePointLocator> octree;
    octree->SetDataSet(polyData);
    octree->BuildLocator();

    for (auto i = 0; i < _N; i++) {
        Vector3d currPos = _positions.col(i);
        vtkNew<vtkIdList> neighbors;
        octree->FindClosestNPoints(2, &(currPos[0]), neighbors);
        neighbors->DeleteId(i);
        _initialNearestNeighbor.push_back(neighbors->GetId(0));
    }
}

//! Function to convert from rotation vectors to normals
void OPSModel::computeNormals(){
    // Assume z-axis of Global Coord Sys is to be rotated
    Quaterniond zaxis(0.0,0.0,0.0,1.0);
    for(auto i=0; i < _N; ++i){
        MapV3d normal( &_normals(0,i), 3, 1 );
        MapV3d rotVec( &_rotationVectors(0,i), 3, 1 );
        Vector3d curr = (Quaterniond(AngleAxisd(rotVec.norm(),
                                                rotVec.normalized()))*zaxis*
                         (Quaterniond(AngleAxisd(rotVec.norm(),
                                                 rotVec.normalized())
                                      ).conjugate())).vec();
        normal << curr(0), curr(1), curr(2);
    }
}

//!Print a VTK file
void OPSModel::printVTKFile(const std::string name){
    vtkNew<vtkCellArray> triangles;
    // Write the finite faces of the triangulation to a VTK file
    for( auto ffi = _dt.finite_faces_begin(); ffi != _dt.finite_faces_end(); ++ffi){
        triangles->InsertNextCell(3);
        for(auto j=2; j >= 0; --j)
            triangles->InsertCellPoint(ffi->vertex(j)->info());
    }

    // Iterate over infinite faces
    Face_circulator fc = _dt.incident_faces(_dt.infinite_vertex()), done(fc);
    if (fc != 0) {
        do{
            triangles->InsertNextCell(3);
            for(auto j=2; j >= 0; --j){
                auto vh = fc->vertex(j);
                auto id = _dt.is_infinite(vh)? 0 : vh->info();
                triangles->InsertCellPoint(id);
            }
        }while(++fc != done);
    }
    //Extract point coordinates for _polyData from x
    vtkNew<vtkDoubleArray> pointCoords;
    pointCoords->SetVoidArray((void*) _positions.data(), 3*_N, 1);
    pointCoords->SetNumberOfComponents(3);

    vtkNew<vtkPoints> points;
    points->SetData( pointCoords );

    //Convert rotation vectors to point normals
    vtkNew<vtkDoubleArray> pointNormals;
    pointNormals->SetName("PointNormals");
    pointNormals->SetVoidArray((void*) _normals.data(), 3*_N, 1);
    pointNormals->SetNumberOfComponents(3);

    // Construct vtkPolyData
    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints( points );
    polyData->GetPointData()->SetNormals(pointNormals);

    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName( name.c_str() );
    writer->SetInputData(polyData);
    writer->Write();
}

//! Calculate average edge length as if the particles were on a mesh
double_t OPSModel::getAverageEdgeLength(){
    double_t avg = 0;
    auto numEdges = 0;
    /*
    for(auto i=0; i < _N; i++){
    Vector3d center = _positions.col(i);
    for(auto j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
    vtkIdType currId = _neighbors[i]->GetId(j);
    Vector3d neighbor = _positions.col(currId);
    avg += (center - neighbor).norm();
    numEdges++;
    }
    }
    */
    avg /= numEdges;
    return avg;
}

//!Compute the OPSBody energy
void OPSModel::compute(){
    // Initialize energies and forces to be zero
    _morseEn = 0.0;
    _normalEn = 0.0;
    _circEn = 0.0;

    // Prepare area constraint variables
    Eigen::Matrix3Xd grad(3,_N);
    grad.setZero(3,_N);
    _value = 0.0;

    computeNormals();
    diffNormalRotVec();

    for(const auto& ring : _first_ring){
        auto i = ring.vertex_id;
        Vector3d p = _normals.col(i);
        Matrix3d M = _diffNormalRV[i];
        std::vector<Vector3d> edgeVecs;

        // Calculate the edge vectors
        for(const auto& j : ring.edges )
            edgeVecs.push_back( _positions.col(j) - _positions.col(i) );

        // Area constraint calculations
        for(const auto& face : ring.faces){
            Vector3d p = edgeVecs[face.first];
            Vector3d q = edgeVecs[face.second];
            double_t S = p.cross(q).norm();
            _value += 0.5*S;
            double_t dAdp = ( q.dot(q)*p - p.dot(q)*q )/(2*S);
            double_t dAdq = ( p.dot(p)*q - p.dot(q)*p )/(2*S);
            grad.col(i) += -1.0*(dAdp + dAdq);
            grad.col(face.first) += dAdp;
            grad.col(face.second) += dAdq;
        }

        for(const auto&j : ring.edges){
            Vector3d q = _normals.col(j);
            Matrix3d N = _diffNormalRV[j];
            Vector3d rn = edgeVecs[j].normalized();
            Vector3d m = p - q;
            Vector3d n = p + q;
            double_t r = edgeVecs[j].norm();
            double_t n_dot_rn = n.dot(rn);

            // Evaluate morse derivatives
            double_t exp_1 = exp( -_a*(r - _re) );
            double_t exp_2 = exp_1*exp_1;
            double_t morseEn =  exp_2 - 2*exp_1;
            Vector3d dMdr = 2*_a*( exp_1 - exp_2 )*rn;

            //Evaluate co-normality derivatives
            double_t Phi_n = m.squaredNorm();
            Vector3d dPhi_nVi = 2*M*m;
            Vector3d dPhi_nVj = -2*N*m;

            //Evaluate co-circularity derivatives
            double_t Phi_c = n_dot_rn*n_dot_rn;
            Vector3d dCdr = (2*n_dot_rn/r)*( n - n_dot_rn*rn );
            Vector3d dPhi_cVi = (2*n_dot_rn)*M*rn;
            Vector3d dPhi_cVj = (2*n_dot_rn)*N*rn;

            // Calculate the total derivatives of energy wrt xi, vi and vj
            Vector3d Dxi = -(dMdr + dCdr/_gamma);
            Vector3d Dvi = (dPhi_nVi + dPhi_cVi )/_gamma;
            Vector3d Dvj = (dPhi_nVj + dPhi_cVj)/_gamma;

            // Update the energies
            _morseEn += morseEn;
            _normalEn += Phi_n/_gamma;
            _circEn += Phi_c/_gamma;
            _f += morseEn + (Phi_n + Phi_c)/_gamma;

            //Update the derivatives
            _posGradient.col(i) += Dxi;
            _posGradient.col(j) -= Dxi; /*!< Neighbor feels opposite force */
            _rotGradient.col(i) += Dvi;
            _rotGradient.col(j) += Dvj;
        }
    }

    double_t areaDiff = _value - _constrainedValue;

    // Add the energy
    _f += 0.5*_K_i*areaDiff*areaDiff - _Lambda_i*areaDiff;

    // Add the gradients
    _posGradient += (_K_i*areaDiff - _Lambda_i)*grad;
}

//! Compute derivative of the normal wrt Rotation Vector
void OPSModel::diffNormalRotVec(){
    for(auto i=0; i < _N; ++i){
        // Read the rotation vector
        Vector3d vi = _rotationVectors.col(i);
        double_t s = sin(0.5*vi.norm()), c = cos(0.5*vi.norm());

        Quaterniond q( AngleAxisd(vi.norm(), vi.normalized()) );
        double_t q0 = q.w(), q1 = q.x(), q2 = q.y(), q3 = q.z();
        Matrix4x3d dpdq;

        dpdq << q2, -q1, q0,
                        q3, -q0, -q1,
                        q0, q3, -q2,
                        q1, q2, q3;
        dpdq = 2 * dpdq;

        Matrix3x4d dqdv;
        dqdv.leftCols(1) = -0.5*s*vi.normalized();
        dqdv.rightCols(3) = (s*Eigen::Matrix3d::Identity() +
                             (0.5*c - s/vi.norm())*vi.normalized()*
                             vi.normalized().transpose())/vi.norm();

        _diffNormalRV[i] = dqdv*dpdq;
    }
}

//! Calculate rotation vector with given point coordinates
void OPSModel::initialRotationVector(RefM3Xd pos, RefM3Xd rotVec){
    //Find unit normal along each point and calculate rotation vector
    // that would map the global z-axis to this unit normal
    for(auto i=0; i < pos.cols(); ++i){
        Vector3d x = pos.col(i);
        x.normalize();
        Vector3d axis, cross_prod;
        double_t angle, cross_prod_norm, p3;

        // Cross-product of x with z-axis.
        cross_prod << -x(1),
                        x(0),
                        0.0;
        cross_prod_norm = cross_prod.norm();
        // Check if x is parallel or anti-parallel to z-axis
        p3 = x(2);
        if( cross_prod_norm < 1e-10){
            axis << 1.0,
                            0.0,
                            0.0; // Arbitrarily choose the x-axis
            angle = (p3 > 0.0)? 0.0 : M_PI;
        }
        else{
            angle = asin( cross_prod_norm );
            angle = (p3 < 0.0)? (M_PI - angle) : angle;
            axis = cross_prod.normalized();
        }
        rotVec.col(i) = angle*axis;
    }
}

//! Minimize Rigid Body motions by applying Kabsch Algorithm
void OPSModel::applyKabschAlgorithm(){
    computeNormals();
    Matrix3Xd pseudoNormal(3,_N);
    pseudoNormal = _positions + _normals;
    Eigen::Affine3d A;
    A = find3DAffineTransform(_positions, _prevX);

    //Apply the transformation to each column in _positions
    for(auto i=0; i < _N; ++i){
        Vector3d tempPos = A.linear()*_positions.col(i) + A.translation();
        Vector3d tempN = A.linear()*pseudoNormal.col(i) + A.translation();
        _positions.col(i) = tempPos;
        _normals.col(i) = (tempN - tempPos).normalized();
    }
    updateRotationVectors();
    //updatePolyData();
}

//! Update Rotation Vectors as per the current normals e.g. after Kabsch update
void OPSModel::updateRotationVectors(){
    for(auto i=0; i < _N; ++i){
        Vector3d x = _normals.col(i);
        Vector3d axis, cross_prod;
        double_t angle, cross_prod_norm, p3;

        // Cross-product of x with z-axis.
        cross_prod << -x(1),
                        x(0),
                        0.0;
        cross_prod_norm = cross_prod.norm();
        // Check if x is parallel or anti-parallel to z-axis
        p3 = x(2);
        if( cross_prod_norm < 1e-10){
            axis << 1.0,
                            0.0,
                            0.0; // Arbitrarily choose the x-axis
            angle = (p3 > 0.0)? 0.0 : M_PI;
        }
        else{
            angle = asin( cross_prod_norm );
            angle = (p3 < 0.0)? (M_PI - angle) : angle;
            axis = cross_prod.normalized();
        }
        _rotationVectors.col(i) = angle*axis;
    }
}

//! Calculate the tangential component of mean-squared displacement
double_t OPSModel::getMeanSquaredDisplacement(){
    _msd_tgt = 0;
    vtkIdType i, nn;
    // We will subtract off the radial displacement.
    for (auto i = 0; i < _N ; ++i) {
        Vector3d xi, xj, diff, xi_diff, xj_diff;
        Vector3d xi0, xj0, xi1, xj1, ni0, nj0;

        nn = _initialNearestNeighbor[i];

        xi0 = _initialPositions.col(i);
        xi1 = _positions.col(i);
        ni0 = xi0.normalized();

        xj0 = _initialPositions.col(nn);
        xj1 = _positions.col(nn);
        nj0 = xj0.normalized();

        xi_diff = (xi1 - xi0);
        xj_diff = (xj1 - xj0);

        xi = xi_diff - (ni0.dot(xi_diff)*ni0);
        xj = xj_diff - (nj0.dot(xj_diff)*nj0);

        diff = xi - xj;
        _msd_tgt += diff.dot(diff);
    }
    _msd_tgt = _msd_tgt / (2*_N);
    return _msd_tgt;
}

//! Calculate mean-squared displacement
//! The tangential MSD calculation has been commented out
std::vector<double_t> OPSModel::getMSD(){
    _msd = 0;
    //_msd_tgt = 0;
    std::vector<double_t> msdAll(2,0.0);
    vtkIdType i, nn;
    // We will subtract off the radial displacement.
    for (auto i = 0; i < _N ; ++i) {
        Vector3d xi, xj, diff, xi_diff, xj_diff;
        Vector3d xi0, xj0, xi1, xj1, ni0, nj0;

        nn = _initialNearestNeighbor[i];

        xi0 = _initialPositions.col(i);
        xi1 = _positions.col(i);
        //ni0 = xi0.normalized();

        xj0 = _initialPositions.col(nn);
        xj1 = _positions.col(nn);
        //nj0 = xj0.normalized();

        xi_diff = (xi1 - xi0);
        xj_diff = (xj1 - xj0);

        //xi = xi_diff - (ni0.dot(xi_diff)*ni0);
        //xj = xj_diff - (nj0.dot(xj_diff)*nj0);

        diff = xi_diff - xj_diff;
        _msd += diff.dot(diff);

        //diff = xi - xj;
        //_msd_tgt += diff.dot(diff);
    }
    _msd = _msd / (2*_N);
    //_msd_tgt = _msd_tgt / _N;
    msdAll[0] = _msd;
    msdAll[1] = 0;
    return msdAll;
}

//! Calculate asphericity
double OPSModel::getAsphericity(){
    double_t asphericity = 0.0;
    double_t R0 = getAverageRadius();
    Eigen::RowVectorXd R(_N);
    R = _positions.colwise().norm();
    asphericity += ((R.array() - R0).square()).sum();
    asphericity /= (_N*R0*R0);
    return asphericity;
}

//! Calculate Volume
double_t OPSModel::getVolume(){
    if(_updateVolume){
        _volume = 0.0;
        /*
    vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
    auto verts = vtkSmartPointer<vtkIdList>::New();
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
    int ida, idb, idc;
    ida = verts->GetId(0);
    idb = verts->GetId(1);
    idc = verts->GetId(2);
    Vector3d a, b, c;
    a = _positions.col(ida);
    b = _positions.col(idb);
    c = _positions.col(idc);

    //Calculate volume
    _volume += (a.dot(b.cross(c)))/6.0;
    }
    _updateVolume = false;
    */
    }
    return _volume;
}

//! Calculate Area
double_t OPSModel::getArea(){
    if(_updateArea){
        _area = 0.0;
        /*
    vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
    auto verts = vtkSmartPointer<vtkIdList>::New();
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
    int ida, idb, idc;
    ida = verts->GetId(0);
    idb = verts->GetId(1);
    idc = verts->GetId(2);
    Vector3d a, b, c;
    a = _positions.col(ida);
    b = _positions.col(idb);
    c = _positions.col(idc);

    //Calculate area
    _area += (b-a).cross(c-a).norm()/2.0;
    }
    _updateArea = false;
    */
    }
    return _area;
}

//! Get average radius
double_t OPSModel::getAverageRadius(){
    if(_updateRadius){
        _radius = 0.0;
        _radius = _positions.colwise().norm().sum()/_N;
        _updateRadius = false;
    }
    return _radius;
}

//! Reset the coordinates of the particles to initial values
//! reset the rotation vectors, polydata and neighbors accordingly
void OPSModel::resetToInitialPositions(){
    _positions = _initialPositions;
    initialRotationVector(_positions, _rotationVectors);
    _posGradient = Matrix3Xd::Zero(3,_N);
    _rotGradient = Matrix3Xd::Zero(3,_N);
    //updatePolyData();
    //updateNeighbors();
}

//! Get Root Mean Squared Angle Deficit
double_t OPSModel::getRMSAngleDeficit(){
    // Prepare a vector to store the angle deficit for each vertex
    Eigen::ArrayXd angleDeficit = VectorXd::Constant( _N, 2*M_PI );
    /*
    // Iterate over all triangles of the polydata
    auto verts = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer< vtkCellArray > faces = _polyData->GetPolys();
    faces->InitTraversal();
    while( faces->GetNextCell( verts ) ){
    // Get all point ids
    vtkIdType id0, id1, id2;
    id0 = verts->GetId(0);
    id1 = verts->GetId(1);
    id2 = verts->GetId(2);
    // Get coordinates of all vertices
    Vector3d v0, v1, v2, e0, e1, e2;
    _polyData->GetPoint( id0, &v0(0) );
    _polyData->GetPoint( id1, &v1(0) );
    _polyData->GetPoint( id2, &v2(0) );
    // Calculate unit vectors along the edges of the triangles
    e0 = (v2 - v1).normalized();
    e1 = (v0 - v2).normalized();
    e2 = (v1 - v0).normalized();
    // Calculate the angles of the triangles
    double_t a0, a1, a2;
    a0 = std::acos( e2.dot( -e1 ) );
    a1 = std::acos( e0.dot( -e2 ) );
    a2 = std::acos( e1.dot( -e0 ) );
    // Subtract the angles from the corresponding vertices' curvature
    angleDeficit(id0) -= a0;
    angleDeficit(id1) -= a1;
    angleDeficit(id2) -= a2;
    }
    // Finally calculate the root mean squared angle deficit
    double_t rmsAngleDeficit = std::sqrt( angleDeficit.square().mean() );
    */
    //return rmsAngleDeficit;
    return 0;
}

void OPSModel::stereoDelaunay(){
    // Copy coordinates
    Matrix3Xd points(3,_N);
    points = _positions;

    // Project points to unit sphere
    points.colwise().normalize();

    // Reset the center of the sphere to origin by translating
    Vector3d center = points.rowwise().mean();
    points = points.colwise() - center;

    // Rotate all points so that the point in 0th column is along z-axis
    Vector3d c = points.col(0);
    double_t cos_t = c(2);
    double_t sin_t = std::sqrt( 1 - cos_t*cos_t );
    Vector3d axis;
    axis << c(1), -c(0), 0.;
    Matrix3d rotMat, axis_cross, outer;
    axis_cross << 0. , -axis(2), axis(1),
                    axis(2), 0., -axis(0),
                    -axis(1), axis(0), 0.;

    outer.noalias() = axis*axis.transpose();

    rotMat = cos_t*Matrix3d::Identity() + sin_t*axis_cross + (1-cos_t)*outer;
    Matrix3Xd rPts(3,_N);
    rPts = rotMat*points; // The points on a sphere rotated

    // Calculate the stereographic projections
    Vector3d p0;
    MapM3Xd l0( &(rPts(0,1)), 3, _N-1 );
    Matrix3Xd l(3,_N-1), proj(3,_N-1);
    p0 << 0,0,-1;
    c = rPts.col(0);
    l = (l0.colwise() - c).colwise().normalized();
    for( auto j=0; j < _N-1; ++j ){
        proj.col(j) = ((p0(2) - l0(2,j))/l(2,j))*l.col(j) + l0.col(j);
    }

    // Insert the projected points in a CGAL vertex_with_info vector
    std::vector< std::pair< Point, unsigned> > verts;
    for( auto j=0; j < _N-1; ++j ){
        verts.push_back(std::make_pair(Point(proj(0,j),proj(1,j)),j+1));
    }

    _dt.clear();
    _dt.insert( verts.begin(), verts.end() );
}

void OPSModel::updateTriangles(){
    stereoDelaunay();
    // *************************** Our data structure *************************//
    _first_ring.clear();
    std::set<std::set<unsigned>> tri, edges;

    // Iterate over all vertices and collect first ring neighbors
    for(auto fvi = _dt.all_vertices_begin(); fvi != _dt.all_vertices_end(); ++fvi){

        vertex_first_ring vfr;
        std::vector<Delaunay::Vertex_handle> rvh;
        auto vid = _dt.is_infinite(fvi)? 0 : fvi->info();
        vfr.vertex_id = vid;
        Delaunay::Edge_circulator ec = _dt.incident_edges(fvi), done(ec);

        // Lambda function to get the vertex id for the edge
        auto getVertexId = [](int a, int b){
            std::set<int> index{0,1,2};
            index.erase(a);
            index.erase(b);
            return *index.begin();
        };

        if( ec != 0){
            do{

                auto fh = ec->first;
                auto edgeIndex = getVertexId(fh->index(fvi),ec->second);
                auto verH = fh->vertex(edgeIndex);
                auto edgeId = _dt.is_infinite(verH)? 0 : verH->info();
                std::set<unsigned> edge{vid,edgeId};
                auto tryInsertEdge = edges.insert(edge);
                if(tryInsertEdge.second){
                    vfr.edges.push_back(edgeId);
                    rvh.push_back(verH);
                }

            }while(++ec != done);
        }

        // Check which edges form a unique face
        auto numEdges = rvh.size();
        for( auto k = 0; k < numEdges; ++k){
            auto next = (k+1) % numEdges;
            // Check if face is formed
            if(_dt.is_face(fvi,rvh[k],rvh[next] )){
                // Check if face is unique
                std::set<unsigned> face{vid,vfr.edges[k],vfr.edges[next]};
                auto tryInsertFace = tri.insert(face);
                if(tryInsertFace.second)
                    vfr.faces.push_back(std::make_pair(k,next));
            }
        }

        _first_ring.push_back(vfr);
    }
    _updateRadius = true;
    _updateVolume = true;
    _updateArea = true;
}

}