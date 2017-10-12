#include "OPSBody.h"

//! Function to update OPSParams
void OPSParams::updateParameter(OPSParams::Parameter p, double_t val){
    switch (p) {
    case OPSParams::D_eV:
        D_e = val;
        break;
    case OPSParams::r_eV:
        r_e = val;
        break;
    case OPSParams::aV:
        a = val;
        break;
    case OPSParams::bV:
        b = val;
        break;
    case OPSParams::gammaV:
        gamma = val;
        break;
    case OPSParams::searchRadiusV:
        searchRadius = val;
        break;
    default:
        std::cout<< "Unknown parameter ignored!" << std::endl;
        break;
    }
}

//! Constructor for BrownOPS
OPSBody::OPSBody(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot, RefM3Xd pG,
                 RefM3Xd rG, OPSParams &p):_f(f),
    _positions(pos.data(),3,n), _rotationVectors(rot.data(),3,n),
    _posGradient(pG.data(),3,n), _rotGradient(rG.data(),3,n),_params(p){

    // Ensure that there is enough memory for n-particles
    assert(n <= pos.cols() && n <= rot.cols() && n <= pG.cols()
           && n <= rG.cols());
    _numPartilces = n;
    _radius = 0.0;
    _volume = 0.0;
    _volConstrained = 0.0;
    _updateRadius = true; /*!< getAverageRadius() will toggle this */
    _updateVolume = true; /*!< getVolume() will toggle this */
    _avgVolConstraintOn = true; /*!< Use average volume constraint by default*/

    _prevX = _positions;

    //Extract point coordinates for _polyData from x
    void *coords = (void*) _positions.data();
    vtkSmartPointer<vtkDoubleArray> pointCoords
            = vtkSmartPointer< vtkDoubleArray >::New();
    pointCoords->SetVoidArray( coords, 3*_numPartilces, 1);
    pointCoords->SetNumberOfComponents(3);

    vtkSmartPointer< vtkPoints > points =
            vtkSmartPointer<vtkPoints>::New();
    points->SetData( pointCoords );

    //Convert rotation vectors to point normals
    _normals = Matrix3Xd::Zero(3,_numPartilces);
    computeNormals();
    void* normalsP = (void*) _normals.data();
    vtkSmartPointer<vtkDoubleArray> pointNormals
            = vtkSmartPointer< vtkDoubleArray >::New();
    pointNormals->SetName("PointNormals");
    pointNormals->SetVoidArray(normalsP, 3*_numPartilces, 1);
    pointNormals->SetNumberOfComponents(3);

    // Construct vtkPolyData
    _polyData = vtkSmartPointer<vtkPolyData>::New();
    _polyData->SetPoints( points );
    _polyData->GetPointData()->SetNormals(pointNormals);

    //Construct the kd-tree
    _kdTree = vtkSmartPointer<vtkKdTree>::New();
    updatePolyDataAndKdTree();

    //Initialize the _neighbors vector
    for(int i=0; i < _numPartilces; i++){
        _neighbors.push_back( vtkSmartPointer<vtkIdList>::New() );
    }
    updateNeighbors();

    //Set the initial nearest neighbor vector
    for (int i = 0; i < _numPartilces; i++) {
        Vector3d currPos = _positions.col(i);
        vtkSmartPointer<vtkIdList> neighbors =
                vtkSmartPointer<vtkIdList>::New();
        _kdTree->FindClosestNPoints(2, &(currPos[0]), neighbors);
        neighbors->DeleteId(i);
        _initialNearestNeighbor.push_back( neighbors->GetId(0) );
    }

}

//! Function to convert from rotation vectors to normals
void OPSBody::computeNormals(){
    // Assume z-axis of Global Coord Sys is to be rotated
    Quaterniond zaxis(0.0,0.0,0.0,1.0);
    for(size_t i=0; i < _numPartilces; ++i){
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

//! Updates _polyData with latest deformed node positions and normals
void OPSBody::updatePolyDataAndKdTree() {
    vtkSmartPointer<vtkPoints> pts =
            vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> unitSphere =
            vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> final;
    vtkSmartPointer<vtkDataSetSurfaceFilter> dssf =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    vtkSmartPointer<vtkIdFilter> idf =
            vtkSmartPointer<vtkIdFilter>::New();
    vtkSmartPointer<vtkDelaunay3D> d3D =
            vtkSmartPointer<vtkDelaunay3D>::New();
    vtkSmartPointer<vtkCellArray> finalCells =
            vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> interim;
    vtkSmartPointer<vtkIdTypeArray> origIds;
    vtkSmartPointer<vtkIdList> pointIds =
            vtkSmartPointer<vtkIdList>::New();

    for(int i=0; i < _numPartilces; i++){
        Vector3d x = _positions.col(i).normalized();
        //We are rounding off to 2 decimals as a hack to aid vtkDelaunay3D
        //Otherwise sometimes we don't get a convex hull
        x[0] = std::round( x[0]*100 )/100;
        x[1] = std::round( x[1]*100 )/100;
        x[2] = std::round( x[2]*100 )/100;
        pts->InsertNextPoint(&(x[0]));
    }
    unitSphere->SetPoints(pts);

    idf->SetIdsArrayName("PointIds");
    idf->PointIdsOn();
    idf->SetInputData(unitSphere);

    //Calculate ideal number of triangles.
    int idealTriCount, pentCount = 12, hexCount;
    hexCount = _numPartilces - pentCount;
    idealTriCount = (6*hexCount + 5*pentCount)/3;

    d3D->SetInputConnection(idf->GetOutputPort());
    dssf->SetInputConnection(d3D->GetOutputPort());
    dssf->Update();
    final = dssf->GetOutput();
    if( final->GetNumberOfPolys() != idealTriCount){
        std::cout<<"The mesh has " << final->GetNumberOfPolys() << " triangles." << std::endl;
        std::cout<< "Bad Delaunay triangulation detected!" <<std::endl;
        /*
        vtkSmartPointer<vtkPolyDataWriter> writer =
                vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(final);
        writer->SetFileName("BadMesh.vtk");
        writer->Write();
        exit(EXIT_FAILURE);*/
    }
    interim = final->GetPolys();
    interim->InitTraversal();
    origIds = vtkIdTypeArray::SafeDownCast(
                final->GetPointData()->GetArray("PointIds"));
    while(interim->GetNextCell(pointIds)){
        int numIds = pointIds->GetNumberOfIds();
        finalCells->InsertNextCell(numIds);
        for(int j=0; j < numIds; j++ ){
            int id = (int)origIds->GetTuple1( pointIds->GetId(j) );
            finalCells->InsertCellPoint(id);
        }
    }
    _polyData->SetPolys(finalCells);
    //Finally update the _kdTree
    _kdTree->BuildLocatorFromPoints( _polyData );

    //Turn on flags to update radius and volume
    _updateRadius = true;
    _updateVolume = true;
    return;
}

//! Store the previous position for Kabsch
void OPSBody::updateDataForKabsch(){
    _prevX = _positions;
}

//! Store the neighbors information for each node
void OPSBody::updateNeighbors(){
    double_t searchR = _params.getSearchRadius();
    for(int i=0; i < _numPartilces; i++){
        Vector3d pos = _positions.col(i);
        _neighbors[i]->Reset();
        _kdTree->FindPointsWithinRadius( searchR, &(pos[0]),
                _neighbors[i] );
        _neighbors[i]->DeleteId(i);
    }
}

//!Print a VTK file
void OPSBody::printVTKFile(const std::string name){
    vtkSmartPointer<vtkPolyDataWriter> writer =
            vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkIdFilter> idf =
            vtkSmartPointer<vtkIdFilter>::New();
    idf->SetInputData(_polyData);
    idf->PointIdsOn();
    idf->SetIdsArrayName("PointIds");
    idf->Update();

    writer->SetFileName( name.c_str() );
    writer->SetInputConnection(idf->GetOutputPort());
    writer->Write();
}

//! Calculate average edge length as if the particles were on a mesh
double_t OPSBody::getAverageEdgeLength(){
    double_t avg = 0;
    int numEdges = 0;
    for(int i=0; i < _numPartilces; i++){
        Vector3d center = _positions.col(i);
        for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
            vtkIdType currId = _neighbors[i]->GetId(j);
            Vector3d neighbor = _positions.col(currId);
            avg += (center - neighbor).norm();
            numEdges++;
        }
    }
    avg /= numEdges;
    return avg;
}

//!Compute the OPSBody energy
void OPSBody::compute(){

    double_t E, re, a, b, G;

    E = _params.getMorseEquilibriumEnergy();
    re = _params.getMorseEquilibriumDistance();
    a = _params.getMorseWellWidth();
    b = _params.getKernelStdDev();
    G = _params.getFVK();

    // Initialize energies and forces to be zero
    _morseEn = 0.0;
    _normalEn = 0.0;
    _circEn = 0.0;

    computeNormals();

    for(int i=0; i < _numPartilces; i++){
        Vector3d vi, xi, p;
        Matrix3d M;

        xi = _positions.col(i);
        vi = _rotationVectors.col(i);
        p = _normals.col(i);
        diffNormalRotVec(vi,M);

        for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
            double_t r, n_dot_rij, exp1, exp2;
            double_t morseEn, Ker, Phi_n, Phi_c;
            Vector3d vj, xj, q, m, n, rij;
            Vector3d dMorseXi, dMorseXj, dKerXi, dKerXj;
            Vector3d dPhi_nVi, dPhi_nVj;
            Vector3d dPhi_cVi, dPhi_cVj;
            Vector3d dPhi_cXi, dPhi_cXj;
            Matrix3d N;

            vtkIdType currId = _neighbors[i]->GetId(j);

            xj = _positions.col(currId);
            vj = _rotationVectors.col(currId);
            q = _normals.col(currId);
            diffNormalRotVec(vj,N);
            rij = xj - xi;
            m = p - q;
            n = p + q;
            r = rij.norm();
            n_dot_rij = n.dot(rij);

            exp1 = exp( -2*a*(r - re) );
            exp2 = exp( -a*(r - re) );

            morseEn = E*( exp1 - 2*exp2 );
            Ker = (E/G)*exp( -r*r/(2*b*b) );
            Phi_n = m.squaredNorm();
            Phi_c = n_dot_rij/r;
            Phi_c *= Phi_c;

            // Evaluate morse derivatives
            dMorseXi = (2*E*a/r)*( exp1 - exp2 )*rij;
            dMorseXj = -dMorseXi;

            // Evaluate kernel derivatives
            dKerXi = (Ker/(b*b))*rij;
            dKerXj = -dKerXi;

            //Evaluate co-normality derivatives
            dPhi_nVi = 2*M*m;
            dPhi_nVj = -2*N*m;

            //Evaluate co-circularity derivatives
            dPhi_cXi = (2*n_dot_rij/(r*r*r*r))*( n_dot_rij*rij -r*r*n );
            dPhi_cXj = -dPhi_cXi;
            dPhi_cVi = (2*n_dot_rij/(r*r))*M*rij;
            dPhi_cVj = (2*n_dot_rij/(r*r))*N*rij;

            Vector3d centerDx, centerDv, neighborDx, neighborDv;

            centerDx = dMorseXi + Ker*dPhi_cXi
                    + dKerXi*(Phi_n + Phi_c );
            centerDv = Ker*(dPhi_nVi + dPhi_cVi );

            neighborDx = dMorseXj + Ker*(dPhi_cXj)
                    + dKerXj*(Phi_n + Phi_c );
            neighborDv = Ker*(dPhi_nVj + dPhi_cVj);

            _posGradient.col(i) += centerDx;
            _posGradient.col(currId) += neighborDx;
            _rotGradient.col(i) += centerDv;
            _rotGradient.col(currId) += neighborDv;

            _morseEn += morseEn;
            _normalEn += Ker*Phi_n;
            _circEn += Ker*Phi_c;
            _f += morseEn + Ker*(Phi_n + Phi_c);
        }
    }
    return;
}

//! Compute derivative of the normal wrt Rotation Vector
void OPSBody::diffNormalRotVec(const RefCV3d &vi, RefM3d diff){
    double_t v0 = vi[0], v1 = vi[1], v2 = vi[2], v = vi.norm();
    double_t s = sin(0.5*v), s_v = s/v, s_v3 = s/(v*v*v);
    double_t c_v2 = 0.5*cos(0.5*v)/(v*v), f = c_v2 - s_v3;

    Quaterniond q( AngleAxisd(v, vi.normalized()) );
    double_t q0 = q.w(), q1 = q.x(), q2 = q.y(), q3 = q.z();
    Matrix3x4d dpdq;

    dpdq << q2, q3, q0, q1,
            -q1, -q0, q3, q2,
            q0, -q1, -q2, q3;
    dpdq = 2 * dpdq;

    Matrix4x3d dqdv;
    Matrix3d V2;
    Vector3d Sc, V2xSc;

    double_t v02 = v0*v0, v12 = v1*v1, v22 = v2*v2;
    double_t v0v1 = v0*v1, v0v2 = v0*v2, v1v2 = v1*v2;

    dqdv.topRows(1) = -0.5*s_v*(vi.transpose());

    V2 << 1, v02, v02,
            1, v12, v12,
            1, v22, v22;
    Sc << s_v,
            -s_v3,
            c_v2;
    V2xSc = V2*Sc;

    dqdv.bottomRows(3) << V2xSc(0), v0v1*f, v0v2*f,
            v0v1*f, V2xSc(1), v1v2*f,
            v0v2*f, v1v2*f, V2xSc(2);

    diff = (dpdq*dqdv).transpose();
}

//! Get Total OPS Energy
double_t OPSBody::getTotalEnergy(){
    return (_morseEn + _normalEn + _circEn);
}

//! Calculate rotation vector with given point coordinates
void OPSBody::initialRotationVector(RefM3Xd pos, RefM3Xd rotVec){
    //Find unit normal along each point and calculate rotation vector
    // that would map the global z-axis to this unit normal
    for(size_t i=0; i < pos.cols(); ++i){
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
void OPSBody::applyKabschAlgorithm(){
    computeNormals();
    Matrix3Xd pseudoNormal(3,_numPartilces);
    pseudoNormal = _positions + _normals;
    Eigen::Affine3d A;
    A = find3DAffineTransform(_positions, _prevX);

    //Apply the transformation to each column in _positions
    for(size_t i=0; i < _numPartilces; ++i){
        Vector3d tempPos = A.linear()*_positions.col(i) + A.translation();
        Vector3d tempN = A.linear()*pseudoNormal.col(i) + A.translation();
        _positions.col(i) = tempPos;
        _normals.col(i) = (tempN - tempPos).normalized();
    }
    updateRotationVectors();
    updatePolyDataAndKdTree();
}

//! Update Rotation Vectors as per the current normals after Kabsch update
void OPSBody::updateRotationVectors(){
    for(size_t i=0; i < _numPartilces; ++i){
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
double_t OPSBody::getMeanSquaredDisplacement(){
    _msd = 0;
    int nn = -1;
    // We will subtract off the radial displacement.
    for (int i = 0; i < _numPartilces; i++) {
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
        _msd += diff.dot(diff);
    }
    _msd = _msd / (2*_numPartilces);
    return _msd;
}

//! Calculate asphericity
double OPSBody::getAsphericity(){
    double asphericity = 0.0;
    double R0 = getAverageRadius();
    Eigen::RowVectorXd R(_numPartilces);
    R = _positions.colwise().norm();
    asphericity += ((R.array() - R0).square()).sum();
    asphericity /= (_numPartilces*R0*R0);
    return asphericity;
}

//! Calculate Volume
double_t OPSBody::getVolume(){
    if(_updateVolume){
        _volume = 0.0;
        vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
        vtkSmartPointer<vtkIdList> verts =
                vtkSmartPointer<vtkIdList>::New();
        cells->InitTraversal();
        while( cells->GetNextCell(verts) ){
            assert(verts->GetNumberOfIds() == 3);
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
    }
    return _volume;
}

//! Get average radius
double_t OPSBody::getAverageRadius(){
    if(_updateRadius){
        _radius = 0.0;
        _radius = _positions.colwise().norm().sum();
        _radius /= _numPartilces;
        _updateRadius = false;
    }
    return _radius;
}

//! Calculate average number of neighbors
double_t OPSBody::getAverageNumberOfNeighbors(){
    double num = 0.0;
    for(int i=0; i < _numPartilces; i++){
        num += _neighbors[i]->GetNumberOfIds();
    }
    num /= _numPartilces;
    return num;
}

//! Calculate the volume constraint and its derivatives
void OPSBody::computeConstraintTerms(double_t &h, RefVXd dh){
    Eigen::Map<Eigen::Matrix3Xd> dhdx(dh.data(),3,_numPartilces);
    if(!_avgVolConstraintOn){
        _volume = 0.0;
        vtkSmartPointer<vtkCellArray> cells = _polyData->GetPolys();
        vtkSmartPointer<vtkIdList> verts =
                vtkSmartPointer<vtkIdList>::New();
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
            _volume += a.dot(b.cross(c))/6.0;

            //Calculate derivative wrt volume
            Vector3d dVa, dVb, dVc;
            dVa <<  b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0];
            dVb << a[2]*c[1]-a[1]*c[2], a[0]*c[2]-a[2]*c[0], a[1]*c[0]-a[0]*c[1];
            dVc << a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0];

            dhdx.col(ida) += 0.1666666666667*dVa;
            dhdx.col(idb) += 0.1666666666667*dVb;
            dhdx.col(idc) += 0.1666666666667*dVc;
        }
    }
    else{
        Eigen::VectorXd R(_numPartilces);
        R = _positions.colwise().norm();
        double_t Ravg = R.sum()/_numPartilces;
        _volume = 4*M_PI*Ravg*Ravg*Ravg/3;
        double_t factor = 4*M_PI*Ravg*Ravg/_numPartilces;
        for(int i=0; i < _numPartilces; ++i){
            dhdx.col(i) = factor*_positions.col(i)/R(i);
        }
    }
    _updateVolume = false;
    h = _volume - _volConstrained;
}

//! Return current value of the volume constraint
double_t OPSBody::getConstraintValue(){
    return (_volume - _volConstrained);
}
