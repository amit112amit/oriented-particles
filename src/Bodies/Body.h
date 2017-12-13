#ifndef BODY_H
#define BODY_H

#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>

namespace OPS{
//! Interface for Body classes
class Body{
public:
    virtual void compute() = 0;
};
}
#endif // BODY_H
