#include "OPSModel.h"

namespace OPS{

//! Constructor for BrownOPS
OPSModel::OPSModel(size_t n, double_t &f, RefVXd x, RefVXd g, RefVXd pX):
        _positions(x.data(),3,n), _rotationVectors(&x(3*n),3,n),
        _posGradient(g.data(),3,n), _rotGradient(&g(3*n),3,n),
        _prevX(pX.data(),3,n), _f(f),_x(x.data(),3*n), _g(g.data(),3*n),
        _pX(pX.data(),3*n){

    // Set number of particles
    _N = n;

    // Initialize internal arrays
    _diffNormalRV = std::vector< Matrix3d,
                    Eigen::aligned_allocator<Matrix3d> >(_N, Matrix3d::Zero());
    //_first_ring.resize(_N);

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
        vtkNew<vtkIdList> neighbors;
        octree->FindClosestNPoints(2, &_positions(0,i), neighbors);
        neighbors->DeleteId(i);
        _initialNearestNeighbor.push_back(neighbors->GetId(0));
    }
    // Brownian body initialization
    std::random_device rd;
    _e2 = std::mt19937(rd());
    _rng = NormD(0.,1.);
    _xi = VectorXd::Zero(3*_N);
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

//! Function to convert from rotation vectors to normals
void OPSModel::computeNormals(const VectorXd &x){
    // Assume z-axis of Global Coord Sys is to be rotated
    Quaterniond zaxis(0.0,0.0,0.0,1.0);
    for(auto i=0; i < _N; ++i){
        MapV3d normal( &_normals(0,i), 3, 1 );
        Vector3d curr = (Quaterniond(AngleAxisd(x.segment<3>(3*(_N+i)).norm(),
                                                x.segment<3>(3*(_N+i)).normalized()))*zaxis*
                         (Quaterniond(AngleAxisd(x.segment<3>(3*(_N+i)).norm(),
                                                 x.segment<3>(3*(_N+i)).normalized())
                                      ).conjugate())).vec();
        normal << curr(0), curr(1), curr(2);
    }
}

//!Print a VTK file
void OPSModel::printVTKFile(const std::string name){
    vtkNew<vtkCellArray> triangles;
    for(const auto& f : _triangles){
        triangles->InsertNextCell(3);
        for(auto j=2; j >= 0; --j)
            triangles->InsertCellPoint(f[j]);
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
    polyData->SetPolys(triangles);

    vtkNew<vtkPolyDataWriter> writer;
    writer->SetFileName( name.c_str() );
    writer->SetInputData(polyData);
    writer->Write();
}

//! Calculate average edge length as if the particles were on a mesh
double_t OPSModel::getAverageEdgeLength(){
    double_t avg = 0;
    for(const auto& e : _edges)
        avg += (_positions.col(e[1]) - _positions.col(e[0])).norm();
    return avg/_edges.size();
}

//!Compute the OPSBody energy
void OPSModel::compute(){
    // Initialize energies and forces to be zero
    _morseEn = 0.0;
    _normalEn = 0.0;
    _circEn = 0.0;

    computeNormals();
    diffNormalRotVec();

    for(const auto& e : _edges){
        // Evaluate morse derivatives
        Vector3d rn = (_positions.col(e[1]) - _positions.col(e[0])).normalized();
        double_t r = (_positions.col(e[1]) - _positions.col(e[0])).norm();
        double_t exp_1 = exp( -_a*(r - _re) );
        double_t exp_2 = exp_1*exp_1;
        double_t morseEn =  exp_2 - 2*exp_1;
        Vector3d dMdr = 2*_a*( exp_1 - exp_2 )*rn;

        //Evaluate co-normality derivatives
        Vector3d m = _normals.col(e[0]) - _normals.col(e[1]);
        double_t Phi_n = m.squaredNorm();
        Matrix3d M = _diffNormalRV[e[0]];
        Vector3d dPhi_nVi = 2*M*m;
        Matrix3d N = _diffNormalRV[e[1]];
        Vector3d dPhi_nVj = -2*N*m;

        //Evaluate co-circularity derivatives
        Vector3d n = _normals.col(e[0]) + _normals.col(e[1]);
        double_t n_dot_rn = n.dot(rn);
        double_t Phi_c = n_dot_rn*n_dot_rn;
        Vector3d dCdr = (2*n_dot_rn/r)*( n - n_dot_rn*rn );
        Vector3d dPhi_cVi = (2*n_dot_rn)*M*rn;
        Vector3d dPhi_cVj = (2*n_dot_rn)*N*rn;

        // Update the energies
        _morseEn += morseEn;
        _normalEn += Phi_n*_gamma_inv;
        _circEn += Phi_c*_gamma_inv;
        _f += morseEn + (Phi_n + Phi_c)*_gamma_inv;

        // Calculate the total derivatives of energy wrt xi, vi and vj
        _posGradient.col(e[0]) -= dMdr + dCdr*_gamma_inv;
        _posGradient.col(e[1]) += dMdr + dCdr*_gamma_inv;
        _rotGradient.col(e[0]) += (dPhi_nVi + dPhi_cVi )*_gamma_inv;
        _rotGradient.col(e[1]) += (dPhi_nVj + dPhi_cVj)*_gamma_inv;
    }

    // Prepare area constraint variables
    Eigen::Matrix3Xd grad(3,_N);
    grad.setZero(3,_N);
    _value = 0.0;
    for(const auto &t : _triangles){
        Vector3d p = _positions.col(t[1]) - _positions.col(t[0]);
        Vector3d q = _positions.col(t[2]) - _positions.col(t[0]);
        double_t S = p.cross(q).norm();
        _value += 0.5*S;
        Vector3d dAdp = ( q.dot(q)*p - p.dot(q)*q )/(2*S);
        Vector3d dAdq = ( p.dot(p)*q - p.dot(q)*p )/(2*S);
        grad.col(t[0]) += -1.0*(dAdp + dAdq);
        grad.col(t[1]) += dAdp;
        grad.col(t[2]) += dAdq;
    }
    double_t areaDiff = _value - _constrainedValue;
    _f += 0.5*_K_i*areaDiff*areaDiff - _Lambda_i*areaDiff;
    _posGradient += (_K_i*areaDiff - _Lambda_i)*grad;

    // Brownian body and viscosity body calculations
    _brownEn = -1.0*_brownCoeff*(_xi.dot(_x - _pX));
    _viscoEn = 0.5*_viscosity*((_x - _pX).dot(_x - _pX));
    _f += _brownEn + _viscoEn;
    _g += -1.0*_brownCoeff*_xi + _viscosity*(_x - _pX);
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
                             (0.5*c - s/vi.norm())*vi*vi.normalized().transpose())/vi.norm();
        _diffNormalRV[i] = dqdv*dpdq;
    }
}

//! Compute derivative of the normal wrt Rotation Vector
void OPSModel::diffNormalRotVec(const VectorXd &x){
    for(auto i=0; i < _N; ++i){
        // Read the rotation vector
        Vector3d vi = x.segment<3>(3*(_N + i));
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
                             (0.5*c - s/vi.norm())*vi*vi.normalized().transpose())/vi.norm();
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

//! Calculate mean-squared displacement
//! The tangential MSD calculation has been commented out
double_t OPSModel::getMSD(){
    _msd = 0;
    vtkIdType nn;
    // We will subtract off the radial displacement.
    for (auto i = 0; i < _N ; ++i) {
        Vector3d diff, xi_diff, xj_diff;
        Vector3d xi0, xj0, xi1, xj1;

        nn = _initialNearestNeighbor[i];

        xi0 = _initialPositions.col(i);
        xi1 = _positions.col(i);

        xj0 = _initialPositions.col(nn);
        xj1 = _positions.col(nn);

        xi_diff = (xi1 - xi0);
        xj_diff = (xj1 - xj0);

        diff = xi_diff - xj_diff;
        _msd += diff.dot(diff);
    }
    _msd /= 2*_N;
    return _msd;
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
        for(const auto & f : _triangles)
            _volume += 0.166666667*(_positions.col(f[0]).
                            dot(_positions.col(f[1]).cross(_positions.col(f[2]))));
        _updateVolume = false;
    }
    return _volume;
}

//! Calculate Area
double_t OPSModel::getArea(){
    if(_updateArea){
        _area = 0.0;
        for(const auto &f : _triangles)
            //Calculate area
            _area += 0.5*(_positions.col(f[1])-_positions.col(f[0])).
                            cross(_positions.col(f[2])-_positions.col(f[0])).norm();
        _updateArea = false;
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
    updateTriangles();
}

//! Get Root Mean Squared Angle Deficit
double_t OPSModel::getRMSAngleDeficit(){
    // Prepare a vector to store the angle deficit for each vertex
    Eigen::ArrayXd angleDeficit = VectorXd::Constant( _N, 2*M_PI );
    for(const auto & f : _triangles){
        // Get coordinates of all vertices
        Vector3d v0 = _positions.col(f[0]);
        Vector3d v1 = _positions.col(f[1]);
        Vector3d v2 = _positions.col(f[2]);

        // Calculate unit vectors along the edges of the triangles
        Vector3d e0 = (v2 - v1).normalized();
        Vector3d e1 = (v0 - v2).normalized();
        Vector3d e2 = (v1 - v0).normalized();

        // Calculate the angles of the triangles
        double_t a0 = std::acos( e2.dot( -e1 ) );
        double_t a1 = std::acos( e0.dot( -e2 ) );
        double_t a2 = std::acos( e1.dot( -e0 ) );

        // Subtract the angles from the corresponding vertices' curvature
        angleDeficit(f[0]) -= a0;
        angleDeficit(f[1]) -= a1;
        angleDeficit(f[2]) -= a2;
    }
    return std::sqrt( angleDeficit.square().mean() );
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
    _edges.clear();
    _triangles.clear();
    // ***************** Our data structure ********************//
    std::set<std::set<unsigned>> tri, edges;

    // Iterate over all vertices and collect first ring neighbors
    for(auto fvi = _dt.all_vertices_begin(); fvi != _dt.all_vertices_end(); ++fvi){

        auto vid = _dt.is_infinite(fvi)? 0 : fvi->info();
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
                    std::array<unsigned,2> ed{vid,edgeId};
                    _edges.push_back(ed);
                }
            }while(++ec != done);
        }

        // Iterate over triangles associated with a vertex
        Delaunay::Face_circulator fc = _dt.incident_faces(fvi), done2(fc);
        if( fc != 0){
            do{
                std::set<unsigned> t;
                std::array<unsigned,3> face;
                for(auto j=0; j < 3; ++j){
                    auto verId = _dt.is_infinite(fc->vertex(j))?
                                            0 : fc->vertex(j)->info();
                    t.insert(verId);
                    face[j] = verId;
                }
                auto tryInsertFace = tri.insert(t);
                if(tryInsertFace.second)
                    _triangles.push_back(face);
            }while(++fc != done2);
        }

    }
    _updateRadius = true;
    _updateVolume = true;
    _updateArea = true;
}

double_t OPSModel::operator()(const VectorXd &x, VectorXd &g){
    MapM3Xd posG(g.data(),3,_N), rotG(&g(3*_N),3,_N);
    g.setZero(6*_N);
    double_t f = 0;
    _morseEn = 0.0;
    _normalEn = 0.0;
    _circEn = 0.0;

    computeNormals(x);
    diffNormalRotVec(x);

    for(const auto& e : _edges){
        // Evaluate morse derivatives
        Vector3d rn = (x.segment<3>(3*e[1]) - x.segment<3>(3*e[0])).normalized();
        double_t r = (x.segment<3>(3*e[1]) - x.segment<3>(3*e[0])).norm();
        double_t exp_1 = exp( -_a*(r - _re) );
        double_t exp_2 = exp_1*exp_1;
        double_t morseEn =  exp_2 - 2*exp_1;
        Vector3d dMdr = 2*_a*( exp_1 - exp_2 )*rn;

        //Evaluate co-normality derivatives
        Vector3d m = _normals.col(e[0]) - _normals.col(e[1]);
        double_t Phi_n = m.squaredNorm();
        Matrix3d M = _diffNormalRV[e[0]];
        Vector3d dPhi_nVi = 2*M*m;
        Matrix3d N = _diffNormalRV[e[1]];
        Vector3d dPhi_nVj = -2*N*m;

        //Evaluate co-circularity derivatives
        Vector3d n = _normals.col(e[0]) + _normals.col(e[1]);
        double_t n_dot_rn = n.dot(rn);
        double_t Phi_c = n_dot_rn*n_dot_rn;
        Vector3d dCdr = (2*n_dot_rn/r)*( n - n_dot_rn*rn );
        Vector3d dPhi_cVi = (2*n_dot_rn)*M*rn;
        Vector3d dPhi_cVj = (2*n_dot_rn)*N*rn;

        // Update the energies
        _morseEn += morseEn;
        _normalEn += Phi_n*_gamma_inv;
        _circEn += Phi_c*_gamma_inv;
        f += morseEn + (Phi_n + Phi_c)*_gamma_inv;

        // Calculate the total derivatives of energy wrt xi, vi and vj
        posG.col(e[0]) -= dMdr + dCdr*_gamma_inv;
        posG.col(e[1]) += dMdr + dCdr*_gamma_inv;
        rotG.col(e[0]) += (dPhi_nVi + dPhi_cVi )*_gamma_inv;
        rotG.col(e[1]) += (dPhi_nVj + dPhi_cVj)*_gamma_inv;
    }

    // Prepare area constraint variables
    Eigen::Matrix3Xd grad(3,_N);
    grad.setZero(3,_N);
    _value = 0.0;
    for(const auto &t : _triangles){
        Vector3d p = x.segment<3>(3*t[1]) - x.segment<3>(3*t[0]);
        Vector3d q = x.segment<3>(3*t[2]) - x.segment<3>(3*t[0]);
        double_t S = p.cross(q).norm();
        _value += 0.5*S;
        Vector3d dAdp = ( q.dot(q)*p - p.dot(q)*q )/(2*S);
        Vector3d dAdq = ( p.dot(p)*q - p.dot(q)*p )/(2*S);
        grad.col(t[0]) += -1.0*(dAdp + dAdq);
        grad.col(t[1]) += dAdp;
        grad.col(t[2]) += dAdq;
    }

    // Area constraint, Brownian body and viscosity body calculations
    double_t areaDiff = _value - _constrainedValue;
    posG += (_K_i*areaDiff - _Lambda_i)*grad;
    _brownEn = -1.0*_brownCoeff*(_xi.dot(x.segment(0,3*_N) - _pX));
    _viscoEn = 0.5*_viscosity*((x.segment(0,3*_N) -_pX).dot(x.segment(0,3*_N)-_pX));
    f += (0.5*_K_i*areaDiff - _Lambda_i)*areaDiff + _brownEn + _viscoEn;
    g.head(3*_N) += -1.0*_brownCoeff*_xi + _viscosity*(x.segment(0,3*_N) - _pX);

    return f;
}

}
