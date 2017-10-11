#if !defined(__VOLCONSTRAINT_H__)
#define __VOLCONSTRAINT_H__

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Body.h"

class VolumeConstraint: public Body{
public:
    typedef Eigen::Matrix3Xd Matrix3Xd;
    typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
    typedef Eigen::Map<Matrix3Xd> Map3Xd;

    VolumeConstraint(size_t N, double_t &f,
                     RefM3Xd pos, RefM3Xd posGrad);
    void compute();
    double_t getAverageVolume();
    double_t getEnergyContribution();
    void setConstrainedVolume(double_t V);
    void updateAugmentedLagrangianCoeffs(double_t L, double_t K);

private:
    size_t _N;
    double_t _avgVol;
    double_t _energyContribution;
    double_t &_f;
    double_t _Lambda_i; /*!< Augmented Lagrangian Multiplier term */
    double_t _k_i; /*!< Augmented Lagrangian penalty stiffness */
    double_t _volConstraint;
    Map3Xd _positions;
    Map3Xd _posGradient;
};

#endif // __VOLCONSTRAINT_H__
