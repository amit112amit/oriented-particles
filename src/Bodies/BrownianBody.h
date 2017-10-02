#if !defined(__BROWNIANBODY_H__)
#define __BROWNIANBODY_H__

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "Body.h"

class BrownianBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;    

    BrownianBody(size_t N, double_t coefficient,
                  double_t &f, RefVXd x, RefVXd g);
    void compute();
    void generateParallelKicks();
    double_t getBrownianEnergy(){return _brownEn;}
    void setCoefficient(double_t C);
private:
    size_t _N;
    double_t _coeff;
    double_t &_f;
    double_t _brownEn;
    MapVXd _x;
    MapVXd _g;
    VectorXd _xi;
    std::mt19937 _e2;
    std::normal_distribution<> _rng;
};

#endif // __BROWNIANBODY_H__

