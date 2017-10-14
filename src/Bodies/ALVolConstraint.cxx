#include "ALVolConstraint.h"

ALVolConstraint::ALVolConstraint(size_t N, double_t &f, RefM3Xd x,
                                 RefM3Xd g):_N(N),_f(f),_xPos(x.data(),3,N),
    _xGrad(g.data(),3,N){
    _volume = 0.0;
    _volConstrained = 0.0;
    _Lambda_i = 1.0;
    _K_i = 1.0;
}

void ALVolConstraint::compute(){
    Eigen::VectorXd R(_N);
    double_t Ravg, factor, volDiff;

    // Calculate avg radius, volume and common factor
    R = _xPos.colwise().norm();
    Ravg = R.sum()/_N;
    _volume = 4.1887902048*Ravg*Ravg*Ravg;
    volDiff = _volume - _volConstrained;

    // Add the energy
    _f += 0.5*_K_i*volDiff*volDiff - _Lambda_i*volDiff;

    factor = 12.5663706144*(_K_i*volDiff - _Lambda_i)*Ravg*Ravg/_N;
    // Add the derivatives
    for(int i=0; i < _N; ++i){
        _xGrad.col(i) += factor*_xPos.col(i)/R(i);
    }
}

void ALVolConstraint::uzawaUpdate(){
    _Lambda_i = _Lambda_i - _K_i*(_volume - _volConstrained);
    _K_i *= 10;
}
