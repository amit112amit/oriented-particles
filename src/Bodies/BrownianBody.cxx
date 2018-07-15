#include "BrownianBody.h"

namespace OPS{
BrownianBody::BrownianBody(size_t N, double_t coeff, double_t &f,
                           RefVXd x, RefVXd g, RefVXd prevX):_N(N),
    _f(f), _x(x.data(),N,1),_g(g.data(),N,1),_prevX(prevX.data(),N,1),
    _coeff(coeff){
    //Initialize the random number generator
    std::random_device rd;
    _e2 = std::mt19937(rd());
    _rng = NormD(0.,1.);
    //Set initial kicks to zero
    _xi = VectorXd::Zero(N);
}

void BrownianBody::generateParallelKicks(){
    for(auto i=0; i < _N; i++){
        _xi(i) = _rng(_e2);
    }
}

void BrownianBody::setCoefficient(double_t C){
    _coeff = C;
}

void BrownianBody::compute(){
    _brownEn = -1.0*_coeff*(_xi.dot(_x - _prevX));
    _f += _brownEn;
    _g += -1.0*_coeff*_xi;
}

void BrownianBody::printVTKFile(const std::string fName){
    void *posPtr = (void*)_x.data();
    auto pts = vtkSmartPointer<vtkPoints>::New();
    auto ptsData = vtkSmartPointer<vtkDoubleArray>::New();
    ptsData->SetVoidArray(posPtr,_N,1);
    ptsData->SetNumberOfComponents(3);
    pts->SetData(ptsData);

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(pts);
    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(fName.c_str());
    writer->SetInputData(poly);
    writer->Write();
}
}
