#include "VolumeConstraint.h"

VolumeConstraint::VolumeConstraint(size_t N, double_t &f, RefM3Xd pos,
                                   RefM3Xd posGrad):_N(N),_f(f),
    _positions(pos.data(),3,N),_posGradient(posGrad.data(),3,N){
    _volConstraint = 0.0;
    _Lambda_i = 1.0;
    _k_i = 1.0;

    // Declare the local variables
    Eigen::RowVectorXd R(_N);
    double_t R0;
    // Calculate individual and average radii
    R = _positions.colwise().norm();
    R0 = R.sum();
    R0 /= _N;

    // Calculate volume and volume difference
    _avgVol = 4.1887902047863905*R0*R0*R0;
}

double_t VolumeConstraint::getAverageVolume(){
    return _avgVol;
}

double_t VolumeConstraint::getEnergyContribution(){
    return _energyContribution;
}

void VolumeConstraint::setConstrainedVolume(double_t V){
    _volConstraint = V;
}

void VolumeConstraint::updateAugmentedLagrangianCoeffs(double_t L, double_t K){
    _Lambda_i = L;
    _k_i = K;
}

void VolumeConstraint::compute(){
    // Declare the local variables
    Eigen::RowVectorXd R(_N);
    double_t R0, factor, volDiff;

    // Calculate individual and average radii
    R = _positions.colwise().norm();
    R0 = R.sum();
    R0 /= _N;

    // Calculate volume and volume difference
    _avgVol = 4.1887902047863905*R0*R0*R0;
    volDiff = _avgVol - _volConstraint;

    // Calculate contribution to the total functional
    _energyContribution = 0.5*_k_i*(volDiff*volDiff) -
            _Lambda_i*(volDiff);
    _f += _energyContribution;

    // Calculate contribution to the derivative wrt x
    factor = 12.5663706144*R0*R0*(_k_i*volDiff -_Lambda_i)/_N;
    for(size_t i=0; i < _N; i++){
        _posGradient.col(i) += (factor/R(i))*_positions.col(i);
    }
}
