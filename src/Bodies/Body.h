#ifndef BODY_H
#define BODY_H

#include <Eigen/Dense>
#include <math.h>
#include <stdio.h>

namespace OPS {
//! Interface for Body classes
class Body {
public:
  virtual void compute() = 0;
};
} // namespace OPS
#endif // BODY_H
