#include "Model.h"

Model::Model(size_t N, double_t &f, RefVXd g):_N(N),_f(f),
    _g(g.data(),N,1){
    std::cout<<"Model DOFs = "<< N << std::endl;
}

void Model::compute(){
    zeroOutData();
    for( auto body : _everyBody){
        body->compute();
    }
}

void Model::addBody(Body *b){
    _everyBody.push_back(b);
}

void Model::zeroOutData(){
    //Zero out the energy and forces;
    _f = 0.0;
    _g.setZero();
}
