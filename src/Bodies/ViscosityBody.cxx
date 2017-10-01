#include "ViscosityBody.h"

//! Constructor
ViscosityBody::ViscosityBody(size_t N, double_t viscosity, double_t &f,
                             const RefCVXd &x, RefVXd g):_N(N),_f(f),
    _x(x),_g(g.data(),1,N),_viscosity(viscosity){
    _prevX = _x;
}

//! Store the current solution as reference for next time step
void ViscosityBody::viscousStep(){
    _prevX = _x;
}

//! Compute the energy and forces
void ViscosityBody::compute(){
    _f += 0.5*_viscosity*((_x-_prevX).dot(_x-_prevX));
    _g += _viscosity*(_x - _prevX);
}
