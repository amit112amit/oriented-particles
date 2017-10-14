#include "ALAreaConstraint.h"

ALAreaConstraint::ALAreaConstraint(size_t N, double_t &f, RefM3Xd x,
                                 RefM3Xd g):_N(N),_f(f),_xPos(x.data(),3,N),
    _xGrad(g.data(),3,N){
    _area = 0.0;
    _areaConstrained = 0.0;
    _Lambda_i = 1.0;
    _K_i = 1.0;
}

void ALAreaConstraint::compute(){
    Eigen::VectorXd R(_N);
    double_t Ravg, factor, areaDiff;

    // Calculate avg radius, volume and common factor
    R = _xPos.colwise().norm();
    Ravg = R.sum()/_N;
    _area = 12.5663706144*Ravg*Ravg;
    areaDiff = _area - _areaConstrained;

    // Add the energy
    _f += 0.5*_K_i*areaDiff*areaDiff - _Lambda_i*areaDiff;

    factor = 25.1327412287*Ravg*(_K_i*areaDiff - _Lambda_i)/_N;
    // Add the derivatives
    for(int i=0; i < _N; ++i){
        _xGrad.col(i) += factor*_xPos.col(i)/R(i);
    }
}

void ALAreaConstraint::uzawaUpdate(){
    _Lambda_i = _Lambda_i - _K_i*(_area - _areaConstrained);
    _K_i *= 10;
}
