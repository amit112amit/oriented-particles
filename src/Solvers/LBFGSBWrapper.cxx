#include "LBFGSBWrapper.h"

#if defined(_WIN32)
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define setulb_ SETULB
#endif
#endif

extern "C" void setulb_(int * n, int *m,
                        double_t * x, double_t * l, double_t * u,
                        int * nbd,
                        double_t * f, double_t * g,
                        double_t * factr, double_t * pgtol,
                        double_t * wa,
                        int * iwa,
                        char * task,
                        int * iprint,
                        char * csave,
                        int * lsave,
                        int * isave,
                        double_t * dsave );

using std::cout;
using std::endl;
using std::setw;
using std::right;
using std::left;
using std::scientific;


//! Constructor
LBFGSBWrapper::LBFGSBWrapper(LBFGSBParams &p, Model &m, double_t &f,
     RefV x, RefV g): _model(m), _f(f), _x(x.data(), x.size(),1),
    _g(g.data(),g.size(),1){
    assert( _x.size() == _g.size());
    _n = _x.size();
    _m = p.getNumHessianCorrections();
    _iprint = p.getPrintCode();
    _maxIterations = p.getMaxIterations();
    _factr = p.getMachineEPSFactor();
    _pgtol = p.getProjectedGradientTolerance();    
    resize(_n);
    _projg = 0.0;
}

void LBFGSBParams::updateProperty(Params p, double_t val){
    switch (p) {
    case m:
        _m = std::round(val);
        break;
    case factr:
        _factr = val;
        break;
    case pgtol:
        _pgtol = val;
        break;
    case iprint:
        _iprint = std::round(val);
        break;
    case maxIterations:
        _maxIterations = std::round(val);
        break;
    default:
        cout<< "Undefined property for LBFGSBParams ignored!" << std::endl;
        break;
    }
}

void LBFGSBWrapper::resize(size_t n){
    _n = n;
    _iwa.resize(3*_n);
    _wa.resize(2*_m*_n+4*_n+12*_m*_m+12*_m);
    _iwa = IntVector_t::Zero(3*_n);
    _wa = Vector_t::Zero(2*_m*_n+4*_n+12*_m*_m+12*_m);
    _nbd.resize(_n);
    _l.resize(_n);
    _u.resize(_n);
    _nbd = IntVector_t::Zero(_n);
    _l = Vector_t::Zero(_n);
    _u = Vector_t::Zero(_n);
}

void LBFGSBWrapper::setBounds(const RefCI nbd, const RefCV l, const RefCV u) {
    if(nbd.size() != _n || l.size() != _n || u.size() != _n ) {
        cout<<"LBFGSBWrapper::setBounds(): input arrays are incorrectly sized."
            << endl;
        return;
    }
    _nbd = nbd;
    _l = l;
    _u = u;
    return;
}


void LBFGSBWrapper::solve() {    
    // set up misc. arrays and data
    char task[60], csave[60];
    for(int i=0; i<60; i++) task[i] = csave[i] = '\0';

    int lsave[4];
    for(int i=0; i<4; i++) lsave[i]=0;

    double_t dsave[29];
    for(int i=0; i<29; i++) dsave[i]=0.0;

    int isave[44];
    for(int i=0; i<44; i++) isave[i]=0;

    cout << "================================================================================"
         << endl << endl
         << "Starting BFGS iterations."
         << endl << endl;

    cout << setw(14) << scientific << right << "|proj g|"
         << setw(14) << scientific << right << "|g|"
         << setw(14) << scientific  << right << "f"
         << setw(14) << right << "iterations"
         << setw(14) << right << "evaluations"
         << endl
         << "--------------------------------------------------------------------------------"
         << endl;

    // We start the iteration by initializing task.

    sprintf(task,"START");

    int iprint=-1;

    _model.compute();

    // ------- the beginning of the loop ----------
    while(true) {

        // This is the call to the L-BFGS-B code.
        setulb_(&_n, &_m, _x.data(), _l.data(), _u.data(), _nbd.data(), &_f,
                _g.data(), &_factr, &_pgtol, _wa.data(),_iwa.data(),
                &(task[0]), &iprint, &(csave[0]),
                &(lsave[0]),&(isave[0]),&(dsave[0]));
        if(strncmp(task,"FG",2)==0) {
            // The minimization routine has returned to request the
            // function f and gradient g values at the current x.

            _model.compute();

            // Go back to the minimization routine.
            continue;
        }

        else if( strncmp(task,"NEW_X",5)==0 ) {
            // stop if maximum number of iterations has been reached
            if ( _maxIterations > 0 && isave[29] > _maxIterations ) {
                break;
            }

            // The minimization routine has returned with a new iterate,
            // and we have opted to continue the iteration.
            if ( _iprint>0 && isave[29]%_iprint == 0 ){
                // 	  char name[100];sprintf(name,"%d",isave[29]);
                cout << setw(14) << scientific << right << dsave[12]
                     << setw(14) << scientific << right << (_g.cwiseAbs()).maxCoeff()
                     << setw(14) << scientific  << right << _f
                     << setw(14) << right << isave[29]
                     << setw(14) << right << isave[33]
                     << endl;
            }

            continue;
        }

        else if( strncmp(task,"CONV",4)==0 ) {
            _model.compute();
            //cout << task << endl;
            std::cout.write(task,49) << std::endl;
            //_model->print("lbfgsbconv");
            break;
        }

        else if( strncmp(task,"ABNORM",6)==0 ) {
            _model.compute();
            std::cout.write(task,49) << std::endl;
            //cout << task << endl;
            //_model->print("abnormal");
            break;
        }

        // If task is neither FG nor NEW_X we terminate execution.
        else {
            //cout << task << endl;
            std::cout.write(task,49) << std::endl;
            break;
        }

        // ---------- the end of the loop -------------
    }
    cout << setw(14) << scientific << right << dsave[12]
         << setw(14) << scientific << right << (_g.cwiseAbs()).maxCoeff()
         << setw(14) << scientific << right << _f
         << setw(14) << right << isave[29]
         << setw(14) << right << isave[33]
         << endl << endl
         << "================================================================================"
         << endl << endl;

    cout.unsetf(std::ios_base::scientific);


    // ---------- the end of solve() -------------
    _projg = dsave[12];
    _iterNo = isave[29];

    _model.compute();
}

//! Solve Augmented Lagrangian Constraints
void LBFGSBWrapper::solveWithALConstraints(AugmentedLagrangian &AL, double_t
                                      conTol, size_t maxIter){
    size_t iter;
    double_t conDiff;

    // Solve the unconstrained minimization problem
    solve();

    conDiff = std::abs(AL.getConstraintValue());
    while( conDiff > conTol && iter < maxIter  ){
        std::cout<< "Augmented Lagrangian Iteration : " << iter << std::endl;

        // Solve the unconstrained minimization problem
        solve();

        // Uzawa-update
        AL.uzawaUpdate();

        // Check termination criteria
        conDiff = std::abs(AL.getConstraintValue());
        iter++;
    }
    std::cout<<"ConDiff = " << conDiff << std::endl;
}