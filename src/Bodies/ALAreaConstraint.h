#if !defined(__ALAREACONSTRAINT_H__)
#define __ALAREACONSTRAINT_H__

#include <Eigen/Dense>
#include "Body.h"

class ALAreaConstraint: public Body{
public:
    typedef Eigen::Matrix3Xd Matrix3Xd;
    typedef Eigen::Map<Matrix3Xd> MapM3Xd;
    typedef Eigen::Ref<Matrix3Xd> RefM3Xd;

    ALAreaConstraint(size_t N, double_t &f, RefM3Xd x,
                                     RefM3Xd g);
    void compute();
    double_t getArea(){return _area;}
    void setConstrainedArea(double_t A){_areaConstrained = A;}
    void setLagrangeCoeff(double_t L){_Lambda_i = L;}
    void setPenaltyCoeff(double_t K){_K_i = K;}
    void uzawaUpdate();

private:
    size_t _N;
    double_t &_f;
    double_t _Lambda_i;
    double_t _K_i;
    double_t _areaConstrained;
    double_t _area;
    MapM3Xd _xPos;
    MapM3Xd _xGrad;
};

#endif //__ALAREACONSTRAINT_H__
