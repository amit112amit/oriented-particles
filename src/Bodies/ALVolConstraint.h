#if !defined(__ALVOLCONSTRAINT_H__)
#define __ALVOLCONSTRAINT_H__

#include <Eigen/Dense>
#include "Body.h"

class ALVolConstraint: public Body{
public:
    typedef Eigen::Matrix3Xd Matrix3Xd;
    typedef Eigen::Map<Matrix3Xd> MapM3Xd;
    typedef Eigen::Ref<Matrix3Xd> RefM3Xd;

    ALVolConstraint(size_t N, double_t &f, RefM3Xd x,
                                     RefM3Xd g);
    void compute();
    double_t getVolume(){return _volume;}
    void setConstrainedVolume(double_t V){_volConstrained = V;}
    void setLagrangeCoeff(double_t L){_Lambda_i = L;}
    void setPenaltyCoeff(double_t K){_K_i = K;}
    void uzawaUpdate();

private:
    size_t _N;
    double_t &_f;
    double_t _Lambda_i;
    double_t _K_i;
    double_t _volConstrained;
    double_t _volume;
    MapM3Xd _xPos;
    MapM3Xd _xGrad;
};

#endif //__ALVOLCONSTRAINT_H__
