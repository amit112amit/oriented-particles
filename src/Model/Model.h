#if !defined(__MODEL_H__)
#define __MODEL_H__

#include <memory>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "Body.h"

namespace OPS{
class Model{
public:
    typedef Eigen::Ref<Eigen::VectorXd> RefVXd;

    Model(size_t N, double_t &f, RefVXd g);
    void addBody(const std::shared_ptr<Body>& b);
    virtual void compute();
    void zeroOutData();
private:
    size_t _N;
    double &_f;
    Eigen::Map<Eigen::VectorXd> _g;
    std::vector<std::shared_ptr<Body>> _everyBody;
};
}
#endif // __MODEL_H__
