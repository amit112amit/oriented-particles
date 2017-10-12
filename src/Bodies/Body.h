#ifndef BODY_H
#define BODY_H

#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>

//! Interface for Body classes
class Body{
public:
    virtual void compute() = 0;
};

//! Interface for Bodies with constraints
class ConstrainedBody{
public:
    typedef Eigen::Ref<Eigen::VectorXd> RefVXd;
    virtual void computeConstraintTerms(double_t &hx, RefVXd dhdx) = 0;
    virtual double_t getConstraintValue() = 0;
};

#endif // BODY_H
