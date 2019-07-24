#include "Pressure.h"

namespace OPS
{

//! Constructor
PressureBody::PressureBody(size_t N, double_t &f, RefM3Xd x, RefM3Xd px,
                           RefM3Xd g, PolyData p)
    : _N(N), _f(f), _x(x.data(), 3, N), _force(g.data(), 3, N), _poly(p),
      _prevX(px.data(), 3, N), _localForce(3, N) {}

//! Compute the energy and forces
void PressureBody::compute()
{
  if (std::abs(_pressure) > 1.0e-8)
  {
    double_t volume = 0.0, prevVol = 0.0;
    _localForce.setZero(3, _N);
    _pvWork = 0;
    auto pts = vtkSmartPointer<vtkIdList>::New();
    auto faces = _poly->GetPolys();
    faces->InitTraversal();
    while (faces->GetNextCell(pts))
    {
      vtkIdType ida, idb, idc;
      Vector3d a, b, c, dVa, dVb, dVc, ap, bp, cp;
      ida = pts->GetId(0);
      idb = pts->GetId(1);
      idc = pts->GetId(2);
      a = _x.col(ida);
      b = _x.col(idb);
      c = _x.col(idc);
      // Calculate volume
      volume += (a.dot(b.cross(c))) / 6.0;
      // Assume that connectivity remains unchanged
      ap = _prevX.col(ida);
      bp = _prevX.col(idb);
      cp = _prevX.col(idc);
      prevVol += (ap.dot(bp.cross(cp))) / 6.0;

      // Calculate derivative wrt volume
      dVa << b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2],
          b[0] * c[1] - b[1] * c[0];
      dVb << a[2] * c[1] - a[1] * c[2], a[0] * c[2] - a[2] * c[0],
          a[1] * c[0] - a[0] * c[1];
      dVc << a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0];

      _localForce.col(ida) += dVa / 6.0;
      _localForce.col(idb) += dVb / 6.0;
      _localForce.col(idc) += dVc / 6.0;
    }
    _localForce *= -_pressure;                /*!< Negative force >*/
    _pvWork = _pressure * (prevVol - volume); /*!< Negative work >*/
    _f += _pvWork;
    _force += _localForce;
  }
}

} // namespace OPS
