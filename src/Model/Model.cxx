#include "Model.h"

Model::Model(size_t N, double_t &f, RefVXd g):_N(N),_f(f),
    _g(g.data(),1,N){}

void Model::compute(){
    //Zero out the energy and forces;
    _f = 0.0;
    _g.setZero();

    for( std::vector<Body*>::iterator i=_everyBody.begin();
         i != _everyBody.end(); ++i){
        (*i)->compute();
    }
}

void Model::addBody(Body *b){
    _everyBody.push_back(b);
}
