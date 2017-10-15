#include "ALConstraint.h"

//! Constructor for Augmented Lagrangian constraint
ALConstraint::ALConstraint(size_t N, double_t &f, RefM3Xd x,
                           RefM3Xd g):_N(N),_f(f),_xPos(x.data(),3,N),
    _xGrad(g.data(),3,N){
    _value = 0.0;
    _constrainedValue = 0.0;
    _Lambda_i = 1.0;
    _K_i = 1.0;
}

//! Uzawa update common to all Augmented Lagrangian constraint subclasses
void ALConstraint::uzawaUpdate(){
    _Lambda_i = _Lambda_i - _K_i*(_value - _constrainedValue);
    _K_i *= 10;
}

//! Calculates average area constraint and derivatives
void AvgAreaConstraint::compute(){
    Eigen::VectorXd R(_N);
    double_t Ravg, factor, areaDiff;

    // Calculate avg radius, volume and common factor
    R = _xPos.colwise().norm();
    Ravg = R.sum()/_N;
    _value = 12.5663706144*Ravg*Ravg;
    areaDiff = _value - _constrainedValue;

    // Add the energy
    _f += 0.5*_K_i*areaDiff*areaDiff - _Lambda_i*areaDiff;

    factor = 25.1327412287*Ravg*(_K_i*areaDiff - _Lambda_i)/_N;
    // Add the derivatives
    for(int i=0; i < _N; ++i){
        _xGrad.col(i) += factor*_xPos.col(i)/R(i);
    }
}

//! Calculates average volume constraint and derivatives
void AvgVolConstraint::compute(){
    Eigen::VectorXd R(_N);
    double_t Ravg, factor, volDiff;

    // Calculate avg radius, volume and common factor
    R = _xPos.colwise().norm();
    Ravg = R.sum()/_N;
    _value = 4.1887902048*Ravg*Ravg*Ravg;
    volDiff = _value - _constrainedValue;


    // Add the energy
    _f += 0.5*_K_i*volDiff*volDiff - _Lambda_i*volDiff;

    factor = 12.5663706144*(_K_i*volDiff - _Lambda_i)*Ravg*Ravg/_N;
    // Add the derivatives
    for(int i=0; i < _N; ++i){
        _xGrad.col(i) += factor*_xPos.col(i)/R(i);
    }
}

//! Calculates exact area constraint and derivatives
void ExactAreaConstraint::compute(){    
    vtkSmartPointer<vtkCellArray> cells = _poly->GetPolys();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    double_t areaDiff, factor;
    Eigen::Matrix3Xd grad(3,_N);
    grad.setZero(3,_N);
    _value = 0.0;
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
        double_t S;
        int ida, idb, idc;
        Eigen::Vector3d a, b, c, p, q, dAdp, dAdq;

        ida = verts->GetId(0);
        idb = verts->GetId(1);
        idc = verts->GetId(2);
        a = _xPos.col(ida);
        b = _xPos.col(idb);
        c = _xPos.col(idc);
        p = b - a;
        q = c - a;
        S = p.cross(q).norm();

        //Calculate area
        _value += S/2.0;

        //Calculate derivatives of area
        dAdp = ( q.dot(q)*p - p.dot(q)*q )/(2*S);
        dAdq = ( p.dot(p)*q - p.dot(q)*p )/(2*S);

        grad.col(ida) += -1.0*(dAdp + dAdq);
        grad.col(idb) += dAdp;
        grad.col(idc) += dAdq;
    }
    areaDiff = _value - _constrainedValue;

    // Add the energy
    _f += 0.5*_K_i*areaDiff*areaDiff - _Lambda_i*areaDiff;

    factor = (_K_i*areaDiff - _Lambda_i);
    // Add the derivatives
    _xGrad += factor*grad;
}

//! Calculates exact volume constraint and derivatives
void ExactVolConstraint::compute(){
    vtkSmartPointer<vtkCellArray> cells = _poly->GetPolys();
    vtkSmartPointer<vtkIdList> verts =
            vtkSmartPointer<vtkIdList>::New();
    double_t volDiff, factor;
    Eigen::Matrix3Xd grad(3,_N);
    grad.setZero(3,_N);
    _value = 0.0;
    cells->InitTraversal();
    while( cells->GetNextCell(verts) ){
        double_t S;
        int ida, idb, idc;
        Eigen::Vector3d a, b, c, dVa, dVb, dVc;

        ida = verts->GetId(0);
        idb = verts->GetId(1);
        idc = verts->GetId(2);
        a = _xPos.col(ida);
        b = _xPos.col(idb);
        c = _xPos.col(idc);

        //Calculate volume
        _value += ( a.dot( b.cross(c) ) )/6.0;

        //Calculate derivative wrt volume
        dVa << b[1]*c[2]-b[2]*c[1],
                b[2]*c[0]-b[0]*c[2],
                b[0]*c[1]-b[1]*c[0];
        dVb << a[2]*c[1]-a[1]*c[2],
                a[0]*c[2]-a[2]*c[0],
                a[1]*c[0]-a[0]*c[1];
        dVc << a[1]*b[2]-a[2]*b[1],
                a[2]*b[0]-a[0]*b[2],
                a[0]*b[1]-a[1]*b[0];

        grad.col(ida) += dVa/6.0;
        grad.col(idb) += dVb/6.0;
        grad.col(idc) += dVc/6.0;
    }
    volDiff = _value - _constrainedValue;

    // Add the energy
    _f += 0.5*_K_i*volDiff*volDiff - _Lambda_i*volDiff;

    factor = (_K_i*volDiff - _Lambda_i);
    // Add the derivatives
    _xGrad += factor*grad;
}
