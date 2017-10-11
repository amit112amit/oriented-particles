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

class BrownianBody: public Body{
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Map<VectorXd> MapVXd;
    typedef Eigen::Ref<VectorXd> RefVXd;    

    BrownianBody(size_t N, double_t coefficient,
                  double_t &f, RefVXd x, RefVXd g, RefVXd prevX);
    void compute();
    void generateParallelKicks();
    double_t getBrownianEnergy(){return _brownEn;}
    void printVTKFile(const std::string fName);
    void setCoefficient(double_t C);

private:
    size_t _N;
    double_t _coeff;
    double_t &_f;
    double_t _brownEn;
    MapVXd _x;
    MapVXd _g;
    MapVXd _prevX;
    VectorXd _xi;    
    std::mt19937 _e2;
    std::normal_distribution<> _rng;
};

#endif // __BROWNIANBODY_H__

