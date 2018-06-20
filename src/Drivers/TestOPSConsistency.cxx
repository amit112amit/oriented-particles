#include <ctime>
#include <stdio.h>
#include <vtkPolyDataReader.h>
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "OPSModel.h"
#include "ViscosityBody.h"
#include "ALConstraint.h"
#include "Pressure.h"

using namespace OPS;

int main(int argc, char* argv[]){

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];
    vtkNew<vtkPolyDataReader> reader;
    reader->SetFileName(inputFileName.c_str());
    reader->Update();
    auto mesh = reader->GetOutput();
    size_t N = mesh->GetNumberOfPoints();
    // ********************************************************//

    // Prepare memory for energy and force
    double_t f;
    Eigen::VectorXd x(6*N), g(6*N), prevX(3*N);
    g.setZero(g.size());
    x.setZero(x.size());

    // Fill arrays
    Eigen::Map<Eigen::Matrix3Xd> pos(x.data(),3,N), rot(&x(3*N),3,N),
            posGrad(g.data(),3,N), rotGrad(&g(3*N),3,N);
    for(auto i = 0; i < N; ++i)
        mesh->GetPoint(i, &pos(0,i));
    prevX = x.head(3*N);
    OPSBody::initialRotationVector(pos, rot);

    // OPS Bodies
    //OPSMesh ops(N,f,pos,rot,posGrad,rotGrad,prevX);
    //ops.updateNeighbors();

    // Create a pressure body
    //PressureBody pBody(N,f,pos,prevPos,posGrad,ops.getPolyData());
    //pBody.setPressure(100.0);

    // Create Brownian and Viscosity bodies
    //Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    //Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    //BrownianBody brown(3*N,1.0,f,thermalX,thermalG,prevX);
    //ViscosityBody visco(3*N,1.0,f,thermalX,thermalG,prevX);

    // Create an Augmented Lagrangian Volume and Area constraint
    //auto poly = ops.getPolyData();
    //auto constraint = ExactAreaConstraint(N,f,pos,posGrad,poly);
    //constraint.setConstraint(1.0);

    // Create Model
    //Model model(6*N,f,g);
    //model.addBody(std::make_shared<OPSMesh>(ops));
    //model.addBody(std::make_shared<ExactAreaConstraint>(constraint));
    //model.addBody(std::make_shared<PressureBody>(pBody));
    //model.addBody(std::make_shared<BrownianBody>(brown));
    //model.addBody(std::make_shared<ViscosityBody>(visco));

    Model model2(6*N,f,g);
    OPSModel opsM(N,f,x,g,prevX);
    opsM.setConstraint(1.0);
    opsM.setLagrangeCoeff(0.0);
    opsM.setPenaltyCoeff(0.0);
    opsM.setViscosity(1.0);
    opsM.setBrownCoeff(1.0);
    model2.addBody(std::make_shared<OPSModel>(opsM));

    // Generate Brownian Kicks
    //brown.generateParallelKicks();
    opsM.generateParallelKicks();

    // Create new x data
    x.setRandom(x.size());
    //x -= 0.2*x;

    // Turn off some bodies
    //brown.setCoefficient(0.0);
    //visco.setViscosity(0.0);
    //constraint.setLagrangeCoeff(0.0);
    //constraint.setPenaltyCoeff(0.0);

// ******************  Consistency Check ********************//

    //Generate log-spaced values
    int hSize = 50;
    double_t start = -12;
    double_t end = -1;
    Eigen::VectorXd hvec(hSize);

    for(auto i=0; i < hSize; ++i){
        double_t currPow = start + i*(end - start)/hSize;
        hvec(i) = std::pow(10,currPow);
    }
    clock_t t = clock();
    // Calculate analytical derivative
    Eigen::VectorXd gAna(g.size());
    //model2.compute();
    opsM(x,g);
    gAna = g;

    // Calculate the error
    for(auto i=0; i < hSize; ++i){
        Eigen::VectorXd gNum(g.size());
        gNum.setZero(g.size());
        double_t h = hvec(i);
        // Calculate the numerical derivative
        for(auto j=0; j < x.size(); ++j){
            double_t fp, fm;
            x(j) += h;
            //model2.compute();
            //fp = f;
            fp = opsM(x,g);
            x(j) += -2*h;
            //model2.compute();
            //fm = f;
            fm = opsM(x,g);
            gNum(j) += (fp-fm)/(2*h);
            x(j) += h;
        }
        std::cout<< hvec(i) <<" , "
                 //<< (gAna - gNum).array().abs().maxCoeff() << std::endl;
                 << (gAna - gNum)<< std::endl;
    }
    std::cout<< "Time elapsed = "
             << ((float)clock() - (float)t)/CLOCKS_PER_SEC
             << " seconds" << std::endl;
    return 1;
}
