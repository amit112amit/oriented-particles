#if !defined(__CONSTRAINTS_H__)
#define __CONSTRAINTS_H__

#include <stdio.h>
#include <math.h>

//! An interface of all constraint types
class Constraint{
public:
    //! Compute the energy and Jacobian terms
    virtual void computeConstraints() = 0;
    virtual double_t getConstraintValue() = 0;
};

#endif //__CONSTRAINTS_H__
