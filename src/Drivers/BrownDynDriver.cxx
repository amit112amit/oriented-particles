#include <stdio.h>
#include <vtkPolyDataReader.h>
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSBody.h"
#include "ViscosityBody.h"

int main(int argc, char* argv[]){

    // Set number of OPS particles
    int N = 5;

    // Prepare memory
    double f;
    Eigen::VectorXd x(6*N), g(6*N);
    x.setRandom(x.size());
    g.setZero(g.size());

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> pos(x.data(),3,N),
            rot(&x(3*N),3,N), posGrad(g.data(),3,N),
            rotGrad(&g(3*N),3,N);
    OPSParams p;
    OPSBody ops(N,f,pos,rot,posGrad,rotGrad,p);

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),1,3*N);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),1,3*N);
    BrownianBody brown(3*N,1.0,f,thermalX,thermalG);
    ViscosityBody visco(3*N,1.0,f,thermalX,thermalG);

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);

    // Create Solver

    return 1;
}
