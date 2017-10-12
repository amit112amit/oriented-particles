#include "Model.h"

Model::Model(size_t N, double_t &f, RefVXd g):_N(N),_f(f),
    _g(g.data(),N,1){
    std::cout<<"Model DOFs = "<< N << std::endl;
}

void Model::compute(){
    zeroOutData();
    for( std::vector<Body*>::iterator i=_everyBody.begin();
         i != _everyBody.end(); ++i){
        (*i)->compute();
    }
    for( std::vector<Constraint*>::iterator i=_constraints.begin();
         i != _constraints.end(); ++i){
        (*i)->computeConstraints();
    }
}

void Model::addConstraint(Constraint *c){
    _constraints.push_back(c);
}

void Model::addBody(Body *b){
    _everyBody.push_back(b);
}

void Model::zeroOutData(){
    //Zero out the energy and forces;
    _f = 0.0;
    _g.setZero();
}
