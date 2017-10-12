#if !defined(__MODEL_H__)
#define __MODEL_H__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "Body.h"
#include "Constraint.h"

class Model{
public:
    typedef Eigen::Ref<Eigen::VectorXd> RefVXd;

    Model(size_t N, double_t &f, RefVXd g);
    void addBody(Body *b);
    void addConstraint(Constraint *c);
    void compute();
    void zeroOutData();
private:
    size_t _N;
    double &_f;
    Eigen::Map<Eigen::VectorXd> _g;
    std::vector<Body*> _everyBody;
    std::vector<Constraint*> _constraints;
};
#endif // __MODEL_H__
