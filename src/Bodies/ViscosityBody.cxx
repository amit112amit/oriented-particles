#include "ViscosityBody.h"

namespace OPS {
//! Constructor
ViscosityBody::ViscosityBody(size_t N, double_t viscosity, double_t &f,
                             RefVXd x, RefVXd g, RefVXd p)
    : _N(N), _f(f), _x(x.data(), N, 1), _g(g.data(), N, 1),
      _prevX(p.data(), N, 1), _viscosity(viscosity) {}

} // namespace OPS
