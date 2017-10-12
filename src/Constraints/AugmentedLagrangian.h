#if !defined(__AUGMENTEDLAGRANGIAN_H__)
#define __AUGMENTEDLAGRANGIAN_H__

#include <stdio.h>
#include <math.h>
#include <Eigen/Dense>
#include "Constraint.h"
#include "Body.h"

class AugmentedLagrangian: public Constraint{

public:
    typedef Eigen::Ref<Eigen::VectorXd> RefVXd;
    AugmentedLagrangian(ConstrainedBody &b, double_t &f, size_t N, RefVXd df);
    void computeConstraints();
    double_t getConstraintValue();
    void setLambdaAndK(double_t L, double_t K);
    void setKFactor(double_t Kfac);
    void uzawaUpdate();

private:
    size_t _N;
    double_t _Lambda_i;
    double_t _K_i;
    double_t _K_factor;
    double_t &_f;
    double_t _hx;
    Eigen::VectorXd _dhdx;
    Eigen::Map<Eigen::VectorXd> _dfdx;
    ConstrainedBody &_body;
};

#endif //__AUGMENTEDLAGRANGIAN_H__
