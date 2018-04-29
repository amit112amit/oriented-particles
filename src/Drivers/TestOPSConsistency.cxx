#include <stdio.h>
#include <vtkPolyDataReader.h>
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"
#include "ALConstraint.h"

using namespace OPS;

int main(int argc, char* argv[]){

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//
    // Set number of OPS particles
    size_t N = mesh->GetNumberOfPoints();

    // Generate Rotation Vectors from input point coordinates
    Eigen::Matrix3Xd coords(3,N);
    for(auto i = 0; i < N; ++i){
        Eigen::Vector3d cp = Eigen::Vector3d::Zero();
        mesh->GetPoint(i, &(cp(0)));
        coords.col(i) = cp;
    }
    Eigen::Matrix3Xd rotVecs(3,N);
    OPSBody::initialRotationVector(coords, rotVecs);

    // Prepare memory for energy and force
    double_t f;
    Eigen::VectorXd x(6*N), g(6*N), prevX(3*N);
    g.setZero(g.size());
    x.setZero(x.size());

    // Fill x with coords and rotVecs and copy coords in prevX
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N), xrot(&(x(3*N)),3,N),
            prevPos(prevX.data(),3,N);
    xpos = coords;
    xrot = rotVecs;
    prevX = x.head(3*N);

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> pos(x.data(),3,N),
            rot(&x(3*N),3,N), posGrad(g.data(),3,N),
            rotGrad(&g(3*N),3,N);    
    OPSMesh ops(N,f,pos,rot,posGrad,rotGrad,prevPos);
    ops.updateNeighbors();

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    BrownianBody brown(3*N,1.0,f,thermalX,thermalG,prevX);
    ViscosityBody visco(3*N,1.0,f,thermalX,thermalG,prevX);

    // Create an Augmented Lagrangian Volume and Area constraint
    vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
    ExactAreaConstraint constraint
            = ExactAreaConstraint(N,f,pos,posGrad,poly);
    constraint.setConstraint(1.0);

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);
    model.addBody(&constraint);

    // Generate Brownian Kicks
    brown.generateParallelKicks();

    // Create new x data
    x.setRandom(x.size());

    // Turn off some bodies
    //brown.setCoefficient(0.0);
    //visco.setViscosity(0.0);
    //constraint.setLagrangeCoeff(0.0);
    //constraint.setPenaltyCoeff(0.0);

// ******************  Consistency Check ********************//

    //Generate log-spaced values
    int hSize = 50;
    double_t start = -12;
    double_t end = -3;
    Eigen::VectorXd hvec(hSize), err(hSize);

    for(auto i=0; i < hSize; ++i){
        double_t currPow = start + i*(end - start)/hSize;
        hvec(i) = std::pow(10,currPow);
    }
    // Calculate analytical derivative
    Eigen::VectorXd gAna(g.size());
    model.compute();

    gAna = g; /*!< Copies the derivative */

    // Calculate the error
    for(auto i=0; i < hSize; ++i){
        Eigen::VectorXd gNum(g.size());
        gNum.setZero(g.size());
        double_t h = hvec(i);
        // Calculate the numerical derivative
        for(auto j=0; j < x.size(); ++j){
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
