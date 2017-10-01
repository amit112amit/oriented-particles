#if !defined(__VISCOSITYBODY_H__)
#define __VISCOSITYBODY_H__

#include <Eigen/Dense>
#include "Body.h"

class ViscosityBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;
    typedef Eigen::Ref<const VectorXd> RefCVXd;

    ViscosityBody(size_t N, double_t viscosity,
                  double_t &f, const RefCVXd &x, RefVXd g);
    void compute();
    void viscousStep();
private:
    size_t _N;
    double_t _viscosity;
    double_t &_f;
    VectorXd _x;
    MapVXd _g;
    VectorXd _prevX;
};

#endif // __VISCOSITYBODY_H__
