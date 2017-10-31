#if !defined(__VISCOSITYBODY_H__)
#define __VISCOSITYBODY_H__

#include <iostream>
#include <Eigen/Dense>
#include "Body.h"

class ViscosityBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;

    ViscosityBody(size_t N, double_t viscosity,
                  double_t &f, RefVXd x, RefVXd g, RefVXd p);
    void compute();
    double_t getViscosityEnergy(){return _viscoEn;}
    VectorXd getPreviousX(){return _prevX;}
    void setViscosity(double_t v);    

private:
    size_t _N;
    double_t _viscosity = 1.0;
    double_t &_f;
    double_t _viscoEn = 0.0;
    MapVXd _x;
    MapVXd _g;
    MapVXd _prevX;
};

#endif // __VISCOSITYBODY_H__
