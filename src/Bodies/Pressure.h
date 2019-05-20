#if !defined(__PRESSURE_H__)
#define __PRESSURE_H__

#include "Body.h"
#include <Eigen/Dense>
#include <iostream>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace OPS {
class PressureBody : public Body {
public:
  typedef Eigen::Vector3d Vector3d;
  typedef Eigen::Matrix3Xd Matrix3Xd;
  typedef Eigen::Map<Matrix3Xd> MapM3Xd;
  typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
  typedef vtkSmartPointer<vtkPolyData> PolyData;

  PressureBody(size_t N, double_t &f, RefM3Xd x, RefM3Xd px, RefM3Xd g,
               PolyData p);
  void compute();
  double_t getPressureWork() { return _pvWork; }
  void getPressureForce(RefM3Xd p) { p = _localForce; }
  void setPressure(double_t p) { _pressure = p; }

private:
  size_t _N;
  double_t _pressure = 0.0;
  double_t &_f;
  double_t _pvWork = 0.0;
  Matrix3Xd _localForce;
  MapM3Xd _x, _prevX;
  MapM3Xd _force;
  PolyData _poly;
};
} // namespace OPS
#endif // __PRESSURE_H__
