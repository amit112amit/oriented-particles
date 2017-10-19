#if !defined(__OPSBODYCUTOFF_H__)
#define __OPSBODYCUTOFF_H__

#include "OPSBody.h"

//! Another version of OPSBody that has a cut-off radius in-built in the potential
//! Note:
//! _r_m > _C + 10*alpha*_N_m.
//! The buffer radius _r_m should satisfy the above. It is driver's responsibility.
//! _N_m is the number of time steps after which to update neighbors.
//! Alpha is the ratio of roughly sqrt(Dt)/r_e where D is Diffusion coefficient
//! t is time step and r_e is equilibrium separation of Morse potential.
class OPSBodyCutOff: public OPSBody{
public:
    OPSBodyCutOff(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot,
                  RefM3Xd posGrad, RefM3Xd rotGrad, OPSParams &p,
                  double_t cutOff, double_t rm);
    void compute();
    void setBufferRadius(double_t b){
        _r_m = b;
    }
    void setCutOff(double_t C){
        _C = C;

    }
    void updateNeighbors();
    void updateRAndD(){
        double_t re = _params.getMorseEquilibriumDistance();
        double_t a = _params.getMorseWellWidth();
        _R = 0.5*(re + (1/a)*log(2) + _C);
        _D = _C - _R;
        _RmD = _R - _D;
        _RpD = _R + _D;
        _PI_2D = M_PI/(2*_D);
        _PI_4D = 0.5*_PI_2D;
        std::cout<< "Updated values of R = "<< _R << " & D = " << _D << std::endl;
    }

private:
    double_t _R;
    double_t _D;
    double_t _C; /*!< Cut Off distance*/
    double_t _RmD;
    double_t _RpD;
    double_t _r_m; /*!< A buffer beyond the cut-off distance */

    // Cache variables to avoid repeated calculations in compute()
    double_t _PI_2D;
    double_t _PI_4D;
};

#endif // __OPSBODYCUTOFF_H__
