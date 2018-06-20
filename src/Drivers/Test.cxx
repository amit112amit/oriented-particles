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

    // Generate Rotation Vectors from input point coordinates
    Eigen::VectorXd x(6*N), g(6*N), prevX(3*N);
    Eigen::Map<Eigen::Matrix3Xd> pos(x.data(),3,N),
            rot(&x(3*N),3,N), posGrad(g.data(),3,N),
            rotGrad(&g(3*N),3,N),prevPos(prevX.data(),3,N);
    for(auto i = 0; i < N; ++i)
        mesh->GetPoint(i, &(pos(0,i)));
    OPSBody::initialRotationVector(pos, rot);
    prevX = x.head(3*N);

    // Prepare memory for energy and force
    double_t f = 0;
    g.setZero(g.size());

    OPSMesh ops(N,f,pos,rot,posGrad,rotGrad,prevPos);
    ops.updateNeighbors();

    // Create an Augmented Lagrangian Volume and Area constraint
    auto poly = ops.getPolyData();
    auto constraint = ExactAreaConstraint(N,f,pos,posGrad,poly);
    constraint.setConstraint(1.0);

    // Create Model
    Model model1(6*N,f,g);
    model1.addBody(std::make_shared<OPSMesh>(ops));
    model1.addBody(std::make_shared<ExactAreaConstraint>(constraint));

    double_t f2 = 0;
    Eigen::VectorXd x2(6*N), g2(6*N), prevX2(3*N);
    x2 = x;
    g2.setZero(g2.size());
    prevX2 = prevX;
    OPSModel opsM(N,f2,x2,g2,prevX2);
    opsM.setConstraint(1.0);
    Model model2(6*N,f2,g2);
    model2.addBody(std::make_shared<OPSModel>(opsM));

    // Create new x data
    x.setRandom(x.size());
    x2 = x;

    clock_t t = clock();
    for(auto i = 0; i < 100000; ++i)
        model1.compute();
    std::cout<< "Old version took " << ((float)clock()-(float)t)/CLOCKS_PER_SEC << std::endl;
    t = clock();
    for(auto i = 0; i < 100000; ++i)
        model2.compute();
    std::cout<< "New version took " << ((float)clock()-(float)t)/CLOCKS_PER_SEC << std::endl;
    for(auto i = 0; i < 100000; ++i)
        opsM(x2,g2);
    std::cout<< "New version took " << ((float)clock()-(float)t)/CLOCKS_PER_SEC << std::endl;
    return 1;
}
