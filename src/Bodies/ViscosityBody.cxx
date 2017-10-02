#include "ViscosityBody.h"

//! Constructor
ViscosityBody::ViscosityBody(size_t N, double_t viscosity, double_t &f,
                             RefVXd x, RefVXd g):_N(N),_f(f),
    _x(x.data(),N,1),_g(g.data(),N,1),_viscosity(viscosity){
    _prevX = _x;
}

//! Store the current solution as reference for next time step
void ViscosityBody::viscousStep(){
    _prevX = _x;
}

//! Compute the energy and forces
void ViscosityBody::compute(){
    _viscoEn = 0.5*_viscosity*((_x-_prevX).dot(_x-_prevX));
    _f += _viscoEn;
    _g += _viscosity*(_x - _prevX);
}

//! Set viscosity value
void ViscosityBody::setViscosity(double_t v){
    _viscosity = v;
}
