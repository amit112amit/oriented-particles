#include "OPSBodyCutOff.h"

//! Constructor of OPSBodyCutOff version
OPSBodyCutOff::OPSBodyCutOff(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot,
                             RefM3Xd posGrad, RefM3Xd rotGrad, OPSParams &p,
                             double_t cutOff, double_t rm):
    OPSBody(n,f,pos,rot,posGrad,rotGrad,p), _C(cutOff), _r_m(rm){
    updateRAndD();
}

//! The main energy computations involving cut-off potential
void OPSBodyCutOff::compute(){
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
    diffNormalRotVec();
    for(int i=0; i < _numPartilces; i++){
        Vector3d vi, xi, p;
        Matrix3d M;

        xi = _positions.col(i);
        vi = _rotationVectors.col(i);
        p = _normals.col(i);
        M = _diffNormalRV[i];

        for(int j=0; j < _neighbors[i]->GetNumberOfIds(); j++){
            Vector3d xj, rij;
            double_t r;

            vtkIdType currId = _neighbors[i]->GetId(j);
            xj = _positions.col(currId);
            rij = xj - xi;
            r = rij.norm();
            if( r >= _RpD)
                continue; /*!< All energies and derivatives will be zero */

            double_t n_dot_rij, exp1, exp2;
            double_t morseEn, Ker, Phi_n, Phi_c, B;
            Vector3d vj, q, m, n;
            Vector3d dMdr, dKdr, dBdr, dPhi_Cdr;
            Vector3d dPhi_nVi, dPhi_nVj;
            Vector3d dPhi_cVi, dPhi_cVj;
            Vector3d centerDx, centerDv, neighborDv;
            Matrix3d N;

            vj = _rotationVectors.col(currId);
            q = _normals.col(currId);
            N = _diffNormalRV[currId];
            m = p - q;
            n = p + q;
            n_dot_rij = n.dot(rij);

            // Evaluate the normality and circularity potentials
            Phi_n = m.squaredNorm();
            Phi_c = n_dot_rij/r;
            Phi_c *= Phi_c;

            //Evaluate co-normality derivatives wrt rotation vectors
            dPhi_nVi = 2*M*m;
            dPhi_nVj = -2*N*m;

            //Evaluate co-circularity derivatives wrt rotation vectors and r
            dPhi_Cdr = (2*n_dot_rij/(r*r*r*r))*( r*r*n - n_dot_rij*rij );
            dPhi_cVi = (2*n_dot_rij/(r*r))*M*rij;
            dPhi_cVj = (2*n_dot_rij/(r*r))*N*rij;

            exp1 = exp( -2*a*(r - re) );
            exp2 = exp( -a*(r - re) );

            // Evaluate morse, kernel and their derivatives
            morseEn = E*( exp1 - 2*exp2 );
            Ker = (E/G)*exp( -r*r/(2*b*b) );
            dMdr = (2*E*a/r)*( exp2 - exp1 )*rij;
            dKdr = (-Ker/(b*b))*rij;
            centerDx = -(dMdr + Ker*dPhi_Cdr + dKdr*(Phi_n + Phi_c ));
            centerDv = Ker*(dPhi_nVi + dPhi_cVi );
            neighborDv = Ker*(dPhi_nVj + dPhi_cVj);
            if( r >= _RmD && r < _RpD){
                // Tersoff switch-off potential
                B = 0.5 - 0.5*sin( _PI_2D*(r - _R) );
                dBdr = -_PI_4D*cos( _PI_2D*(r - _R) )*rij/r;
                dMdr = morseEn*dBdr + B*dMdr;
                dKdr = Ker*dBdr + B*dKdr;
                centerDx = B*centerDx - (morseEn + Ker*(Phi_n + Phi_c))*dBdr;
                centerDv *= B;
                neighborDv *= B;
                morseEn *= B;
                Ker *= B;
            }            

            _posGradient.col(i) += centerDx;
            _posGradient.col(currId) -= centerDx; /*!< neighbor feels opposite force */
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

//! Store the neighbors information for each node
void OPSBodyCutOff::updateNeighbors(){
    for(int i=0; i < _numPartilces; i++){
        Vector3d pos = _positions.col(i);
        _neighbors[i]->Reset();
        // Use buffer radius instead of just cut-off
        _octree->FindPointsWithinRadius( _r_m, &(pos[0]),
                _neighbors[i] );
        _neighbors[i]->DeleteId(i);
    }
    // To avoid double counting we should remove the repeated
    // neighbors
    for(int i=0; i < _numPartilces; ++i){
        for(int j=0; j < _neighbors[i]->GetNumberOfIds(); ++j){
            vtkIdType currId = _neighbors[i]->GetId(j);
            _neighbors[currId]->DeleteId(i);
        }
    }
}
