#include <stdio.h>
#include <vtkPolyDataReader.h>
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSBody.h"
#include "ViscosityBody.h"
#include "VolumeConstraint.h"

int main(int argc, char* argv[]){

    // Set number of OPS particles
    int N = 5;

    // Prepare memory
    double_t f;
    Eigen::VectorXd x(6*N+1), g(6*N+1), prevX(3*N);
    x.setRandom(x.size());
    g.setZero(g.size());
    prevX.setRandom(prevX.size());

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> pos(x.data(),3,N),
            rot(&x(3*N),3,N), posGrad(g.data(),3,N),
            rotGrad(&g(3*N),3,N);
    OPSParams p;
    OPSBody ops(N,f,pos,rot,posGrad,rotGrad,p);

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    BrownianBody brown(3*N,1.0,f,thermalX,thermalG,prevX);
    ViscosityBody visco(3*N,1.0,f,thermalX,thermalG,prevX);

    // Create Volume constraint body
    VolumeConstraint volC(N,f,pos,posGrad);

    // Create Model
    Model model(6*N+1,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);
    model.addBody(&volC);

    // Generate Brownian Kicks
    brown.generateParallelKicks();

    // Create new x data
    x.setRandom(x.size());

    // Turn off some bodies
    //p.updateParameter(OPSParams::D_eV,0.0);
    //brown.setCoefficient(0.0);
    //visco.setViscosity(0.0);
    //volC.updateAugmentedLagrangianCoeffs(0.0,0.0);

// ******************  Consistency Check ********************//

    //Generate log-spaced values
    int hSize = 50;
    double_t start = -12;
    double_t end = -3;
    Eigen::VectorXd hvec(hSize), err(hSize);

    for(int i=0; i < hSize; ++i){
        double_t currPow = start + i*(end - start)/hSize;
        hvec(i) = std::pow(10,currPow);
    }
    // Calculate analytical derivative
    Eigen::VectorXd gAna(g.size());
    model.compute();
    std::cout <<"OPSBody energy = " << ops.getTotalEnergy() << std::endl;
    std::cout <<"BrownBody energy = " << brown.getBrownianEnergy() << std::endl;
    std::cout <<"ViscoBody energy = " << visco.getViscosityEnergy() << std::endl;
    std::cout <<"VolConstr energy = " << volC.getEnergyContribution() << std::endl;
    gAna = g; /*!< Copies the derivative */

    // Calculate the error
    for(int i=0; i < hSize; ++i){
        Eigen::VectorXd gNum(g.size());
        gNum.setZero(g.size());
        double_t h = hvec(i);
        // Calculate the numerical derivative
        for(int j=0; j < x.size(); ++j){
            double_t fp, fm;
            x(j) += h;
            model.compute();
            fp = f;
            x(j) += -2*h;
            model.compute();
            fm = f;
            gNum(j) += (fp-fm)/(2*h);
            x(j) += h;
        }
        err(i) = (gAna - gNum).array().abs().maxCoeff();
        std::cout<< hvec(i) <<" , " << err(i) << std::endl;
    } 
    return 1;
}
