#if !defined(__ALCONSTRAINT_H__)
#define __ALCONSTRAINT_H__

#include <Eigen/Dense>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include "Body.h"

//! Parent class for all Augmented Lagrangian constraints
class ALConstraint: public Body{
public:
    typedef Eigen::Matrix3Xd Matrix3Xd;
    typedef Eigen::Map<Matrix3Xd> MapM3Xd;
    typedef Eigen::Ref<Matrix3Xd> RefM3Xd;

    ALConstraint(size_t N, double_t &f, RefM3Xd x, RefM3Xd g);
    virtual void compute() = 0;
    double_t getConstraintValue(){return (_value - _constrainedValue);}
    void setConstraint(double_t V){_constrainedValue = V;}
    void setLagrangeCoeff(double_t L){_Lambda_i = L;}
    void setPenaltyCoeff(double_t K){_K_i = K;}
    void uzawaUpdate();

protected:
    size_t _N;
    double_t &_f;
    double_t _Lambda_i;
    double_t _K_i;
    double_t _constrainedValue;
    double_t _value;
    MapM3Xd _xPos;
    MapM3Xd _xGrad;
};

//! Implements an average area Augmented Lagrangian constraint
class AvgAreaConstraint: public ALConstraint{
public:
    AvgAreaConstraint(size_t N, double_t &f, RefM3Xd x,
                      RefM3Xd g):ALConstraint(N,f,x,g){}
    void compute();
};

//! Implements an average volume Augmented Lagrangian constraint
class AvgVolConstraint: public ALConstraint{
public:
    AvgVolConstraint(size_t N, double_t &f, RefM3Xd x,
                      RefM3Xd g):ALConstraint(N,f,x,g){}
    void compute();
};

//! Implements an exact area Augmented Lagrangian constraint
class ExactAreaConstraint: public ALConstraint{
public:
    ExactAreaConstraint(size_t N, double_t &f, RefM3Xd x,
                      RefM3Xd g, vtkSmartPointer<vtkPolyData> p):
        ALConstraint(N,f,x,g){
        _poly = p;
    }
    void compute();
private:
    vtkSmartPointer<vtkPolyData> _poly;
};

//! Implements an exact volume Augmented Lagrangian constraint
class ExactVolConstraint: public ALConstraint{
public:
    ExactVolConstraint(size_t N, double_t &f, RefM3Xd x,
                      RefM3Xd g, vtkSmartPointer<vtkPolyData> p):
        ALConstraint(N,f,x,g){
        _poly = p;
    }
    void compute();
private:
    vtkSmartPointer<vtkPolyData> _poly;
};
#endif //__ALCONSTRAINT_H__
