#include "AugmentedLagrangian.h"

//! Constructor takes a ConstrainedBody and location to update energy and
//! Jacobian
AugmentedLagrangian::AugmentedLagrangian(ConstrainedBody &b, double_t &f,
                                         size_t N, RefVXd df):_body(b), _f(f),
    _N(N),_dfdx(df.data(),N,1){
    _Lambda_i = 1.0;
    _K_i = 100.0;
    _K_factor = 10.0;
    _hx = 0.0;
    _dhdx.setZero(N);
}

//! Update the factor by which to update the penalty coefficient
void AugmentedLagrangian::setKFactor(double_t Kfac){
    _K_factor = Kfac;
}

//! Set the values of Lambda and K to non-default
void AugmentedLagrangian::setLambdaAndK(double_t L, double_t K){
    _Lambda_i = L;
    _K_i = K;
}

//! Update Lambda and K as per Uzawa algorithm
void AugmentedLagrangian::uzawaUpdate(){
    _Lambda_i = _Lambda_i - _K_i*_hx;
    _K_i = _K_factor*_K_i;
}

//! Compute the energy and gradient contributions of constraint
void AugmentedLagrangian::computeConstraints(){
    _body.computeConstraintTerms(_hx,_dhdx);
    _f += 0.5*_K_i*_hx*_hx - _Lambda_i*_hx;
    _dfdx += (_K_i*_hx - _Lambda_i)*_dhdx;
}

//! Get the termination values from the constrained body
double_t AugmentedLagrangian::getConstraintValue(){
    return _body.getConstraintValue();
}
