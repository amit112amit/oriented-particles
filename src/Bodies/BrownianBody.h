#if !defined(__BROWNIANBODY_H__)
#define __BROWNIANBODY_H__

#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "Body.h"
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

namespace OPS{
class BrownianBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;
    typedef std::normal_distribution<double_t> NormD;

    BrownianBody(size_t N, double_t coefficient,
                 double_t &f, RefVXd x, RefVXd g, RefVXd prevX);
    inline void compute(){
        _brownEn = -1.0*_coeff*(_xi.dot(_x - _prevX));
        _f += _brownEn;
        _g += -1.0*_coeff*_xi;
    }
    inline void generateParallelKicks(){
        for(auto i=0; i < _N; i++){
        _xi(i) = _rng(_e2);
        }
    }
    inline std::mt19937 getRandomEngine(){return _e2;}
    inline NormD getRandomGenerator(){return _rng;}
    inline double_t getBrownianEnergy(){return _brownEn;}
    void printVTKFile(const std::string fName);
    inline void setCoefficient(double_t C){_coeff = C;}
    inline void setRandomEngine( std::mt19937 e ){_e2 = e;}
    inline void setRandomGenerator( NormD r ){_rng = r;}

private:
    size_t _N;
    double_t _coeff = 1.0;
    double_t &_f;
    double_t _brownEn = 0.0;
    MapVXd _x;
    MapVXd _g;
    MapVXd _prevX;
    VectorXd _xi;
    std::mt19937 _e2;
    NormD _rng;
};
}
#endif // __BROWNIANBODY_H__

