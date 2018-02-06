#include "PeriodicPlateBody.h"

namespace OPS{

    // Constructor
    PeriodicPlateBody::PeriodicPlateBody(size_t N, double_t &f, Matrix3Xd &p,
	    Matrix3Xd &g, Matrix3Xd &r, Matrix3Xd &rg, double_t L, double_t W):
	_N(N), _f(f), _positions(p.data(),3,N), _posGradient(g.data(),3,N),
	_rotationVectors(r.data(),3,N), _rotGradient(rg.data(),3,N), _L(L),
	_W(W), _polyData( vtkSmartPointer<vtkPolyData>::New() ),
	_kdtree( vtkSmartPointer<vtkKdTreePointLocator>::New() ){
	    // Ensure that there is enough memory for n-particles
	    assert(N <= pos.cols() && N <= rot.cols() && N <= pG.cols()
		    && N <= rG.cols());

	    // Set number of particles
	    _N = N;

	    // Initialize internal arrays
	    _prevX = _positions;
	    _diffNormalRV = std::vector< Matrix3d >(_N,Matrix3d::Zero());

	    //Extract point coordinates for _polyData from x
	    void *coords = (void*) _positions.data();
	    auto pointCoords = vtkSmartPointer< vtkDoubleArray >::New();
	    pointCoords->SetVoidArray( coords, 3*_N, 1);
	    pointCoords->SetNumberOfComponents(3);

	    auto points = vtkSmartPointer<vtkPoints>::New();
	    points->SetData( pointCoords );

	    //Convert rotation vectors to point normals
	    _normals = Matrix3Xd::Zero(3,_N);
	    computeNormals();
	    void* normalsP = (void*) _normals.data();
	    auto pointNormals = vtkSmartPointer< vtkDoubleArray >::New();
	    pointNormals->SetName("PointNormals");
	    pointNormals->SetVoidArray(normalsP, 3*_N, 1);
	    pointNormals->SetNumberOfComponents(3);

	    // Construct vtkPolyData
	    _polyData->SetPoints( points );
	    _polyData->GetPointData()->SetNormals(pointNormals);

	    //Construct the kd-tree
	    updatePeriodicKdTree();

	    //Initialize the _neighbors vector
	    for(auto i=0; i < _N; i++){
		_neighbors.push_back( vtkSmartPointer<vtkIdList>::New() );
	    }
	    updateNeighbors();

	    //Set the initial nearest neighbor vector
	    for (auto i = 0; i < _N; i++) {
		Vector3d currPos = _positions.col(i);
		auto neighbors = vtkSmartPointer<vtkIdList>::New();
		_kdtree->FindClosestNPoints(2, &(currPos[0]), neighbors);
		neighbors->DeleteId(i);
		_initialNearestNeighbor.push_back( neighbors->GetId(0) );
	    }
	} 

    // The compute function
    void PeriodicPlateBody::compute(){

	// Initialize energies and forces to be zero
	_morseEn = 0.0;
	_normalEn = 0.0;
	_circEn = 0.0;

	computeNormals();
	diffNormalRotVec();
	for(auto i=0; i < _N; i++){
	    Vector3d vi, xi, p;
	    Matrix3d M;

	    xi = _positions.col(i);
	    vi = _rotationVectors.col(i);
	    p = _normals.col(i);
	    M = _diffNormalRV[i];

	    for(auto j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
		double_t r, n_dot_rij, exp_2, exp_1;
		double_t morseEn, Ker, Phi_n, Phi_c;
		Vector3d vj, xj, q, m, n, rij;
		Vector3d dMorseXi, dMorseXj, dKerXi, dKerXj;
		Vector3d dPhi_nVi, dPhi_nVj;
		Vector3d dPhi_cVi, dPhi_cVj;
		Vector3d dPhi_cXi, dPhi_cXj;
		Matrix3d N;

		vtkIdType currId = _neighbors[i]->GetId(j);
		vtkIdType mappedId = periodicMap[currId];

		xj = _positions.col(currId);
		vj = _rotationVectors.col(mappedId);
		q = _normals.col(mappedId);
		N = _diffNormalRV[mappedId];
		rij = xj - xi;
		m = p - q;
		n = p + q;
		r = rij.norm();
		n_dot_rij = n.dot(rij);

		exp_1 = exp( -_a*(r - _re) );
		exp_2 = exp_1*exp_1;

		//morseEn = _De*( exp_2 - 2*exp_1 );
		morseEn =  exp_2 - 2*exp_1;
		//Ker = (_De/_gamma)*exp( -r*r/(2*_b*_b) );
		Phi_n = m.squaredNorm();
		Phi_c = n_dot_rij/r;
		Phi_c *= Phi_c;

		// Evaluate morse derivatives
		//dMorseXi = (2*_De*_a/r)*( exp_2 - exp_1 )*rij;
		dMorseXi = (2*_a/r)*( exp_2 - exp_1 )*rij;
		dMorseXj = -dMorseXi;

		// Evaluate kernel derivatives
		//dKerXi = (Ker/(_b*_b))*rij;
		//dKerXj = -dKerXi;

		//Evaluate co-normality derivatives
		dPhi_nVi = 2*M*m;
		dPhi_nVj = -2*N*m;

		//Evaluate co-circularity derivatives
		dPhi_cXi = (2*n_dot_rij/(r*r*r*r))*( n_dot_rij*rij -r*r*n );
		dPhi_cXj = -dPhi_cXi;
		dPhi_cVi = (2*n_dot_rij/(r*r))*M*rij;
		dPhi_cVj = (2*n_dot_rij/(r*r))*N*rij;

		Vector3d centerDx, centerDv, neighborDx, neighborDv;
		/*
		   centerDx = dMorseXi + Ker*dPhi_cXi
		   + dKerXi*(Phi_n + _circCoeff*Phi_c );
		   centerDv = Ker*(dPhi_nVi + _circCoeff*dPhi_cVi );

		   neighborDx = dMorseXj + Ker*(dPhi_cXj)
		   + dKerXj*(Phi_n + _circCoeff*Phi_c );
		   neighborDv = Ker*(dPhi_nVj + _circCoeff*dPhi_cVj);
		   */
		centerDx = dMorseXi + dPhi_cXi/_gamma;
		centerDv = (dPhi_nVi + _circCoeff*dPhi_cVi )/_gamma;

		neighborDx = dMorseXj + (dPhi_cXj)/_gamma;
		neighborDv = (dPhi_nVj + _circCoeff*dPhi_cVj)/_gamma;

		_posGradient.col(i) += centerDx;
		_posGradient.col(mappedId) += neighborDx;
		_rotGradient.col(i) += centerDv;
		_rotGradient.col(mappedId) += neighborDv;

		_morseEn += morseEn;
		_normalEn += Phi_n/_gamma;
		_circEn += _circCoeff*Phi_c/_gamma;

		_f += morseEn + (Phi_n + _circCoeff*Phi_c)/_gamma;
	    }
	}
	return;
    }

    // The crucial update Kd-Tree function
    PeriodicPlateBody::updatePeriodicKdTree(){

    }
}
