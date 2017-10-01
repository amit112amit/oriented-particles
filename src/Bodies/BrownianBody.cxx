#include "BrownianBody.h"

BrownianBody::BrownianBody(size_t N, double_t coeff, double_t &f,
                           const RefCVXd &x, RefVXd g):_N(N), _f(f),
    _x(x),_g(g.data(),1,N),_coeff(1.41421356237*coeff){
    //Initialize the random number generator
    std::random_device rd;
    _e2 = std::mt19937(rd());
    _rng = std::normal_distribution<>(0,1);
    //Set initial kicks to zero
    _xi = VectorXd::Zero(N);
}

void BrownianBody::generateParallelKicks(){
    for(int i=0; i < _N; i++){
        _xi(i) = _rng(_e2);
    }
}

void BrownianBody::setCoefficient(double_t C){
    _coeff = 1.41421356237*C;
}

void BrownianBody::compute(){
    _f -= _coeff*(_xi.dot(_x));
    _g -= _coeff*_xi;
}
