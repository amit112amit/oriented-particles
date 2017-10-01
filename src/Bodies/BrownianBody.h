#if !defined(__BROWNIANBODY_H__)
#define __BROWNIANBODY_H__

#include <random>
#include <Eigen/Dense>
#include "Body.h"

class BrownianBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;
    typedef Eigen::Ref<const VectorXd> RefCVXd;

    BrownianBody(size_t N, double_t coefficient,
                  double_t &f, const RefCVXd &x, RefVXd g);
    void compute();
    void generateParallelKicks();
    void setCoefficient(double_t C);
private:
    size_t _N;
    double_t _coeff;
    double_t &_f;
    VectorXd _x;
    MapVXd _g;
    VectorXd _xi;
    std::mt19937 _e2;
    std::normal_distribution<> _rng;
};

#endif // __BROWNIANBODY_H__

