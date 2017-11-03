// Reference:
//  1. Nonlinear quasi-Newton solver using L-BFGS-B code of Jorge Nocedal.
//  http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
//  2. https://github.com/wsklug/voom -> src/Solvers/Lbfgsb.h
//
/////////////////////////////////////////////////////////////////////////

#if !defined(__LBFGSBWRAPPER_H__)
#define __LBFGSBWRAPPER_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ios>
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "Model.h"

#ifdef _NO_PRINTING_
#define PRINT(X) do {} while(0)
#define PRINT_CHAR_ARR(X,N) do {} while(0)
#define UNSET(X) do {} while(0)
#else
#define PRINT(X) do { std::cout << X; } while(0)
#define PRINT_CHAR_ARR(X,N) do { std::cout.write(X,N) << std::endl; } while(0)
#define UNSET(X) do { std::cout.unsetf(X); } while(0)
#endif

//! l-BFGS-b Solver Parameter Class
class LBFGSBParams{
public:
    //! Default Constructor
    LBFGSBParams():_m(5),_iprint(1000),_maxIterations(100000),
    _factr(10.0),_pgtol(1e-8){}

    //! Constructor with parameters
    LBFGSBParams(size_t m = 5, size_t i = 0, size_t mx = 1000,
                 double_t f = 10, double_t p = 1.0e-5):_m(m),
        _iprint(i), _maxIterations(mx),_factr(f), _pgtol(p){}

    enum Params{m,factr,pgtol,iprint,maxIterations};

    void updateProperty(Params p, double_t val);    
    size_t getNumHessianCorrections(){return _m;}
    size_t getPrintCode(){return _iprint;}
    size_t getMaxIterations(){return _maxIterations;}
    double_t getMachineEPSFactor(){return _factr;}
    double_t getProjectedGradientTolerance(){return _pgtol;}

private:
    int _n, _m, _iprint, _maxIterations;
    double _factr, _pgtol;
};

//! The l-BFGS-b Wrapper class
class LBFGSBWrapper{
public:
    typedef Eigen::VectorXd Vector_t;
    typedef Eigen::Ref<Vector_t> RefV;
    typedef Eigen::Ref<const Vector_t> RefCV;
    typedef Eigen::Map<Vector_t> MapV;
    typedef Eigen::VectorXi IntVector_t;
    typedef Eigen::Ref<IntVector_t> RefI;
    typedef Eigen::Ref<const IntVector_t> RefCI;

    LBFGSBWrapper(LBFGSBParams &p, Model &s, double_t &f, RefV x, RefV g);
    void setBounds(const RefCI nbd, const RefCV l, const RefCV u);
    void setNumHessianCorrections(size_t m){ _m=m; }
    void setPrintCode(size_t p){ _iprint=p; }
    void setMaxIterations(size_t i){ _maxIterations=i; }
    void setMachineEPSFactor(double_t f){ _factr=f; }
    void setProjectedGradientTolerance(double_t p){ _pgtol=p; }
    void solve();
    void turnOffLogging(){ _loggingOn = false;}

private:
    void resize(size_t n);
    bool _loggingOn = true;
    double_t &_f;
    MapV _x;
    MapV _g;
    Model &_model;    
    Vector_t _l;
    Vector_t _u;
    Vector_t _wa;
    IntVector_t _nbd;
    IntVector_t _iwa;

    int _n, _m, _iprint, _maxIterations, _iterNo;
    double _factr, _pgtol;
    double _projg;

};
#endif // __LBFGSBWRAPPER_H__
