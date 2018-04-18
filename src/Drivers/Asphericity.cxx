#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"

using namespace OPS;

int main(int argc, char* argv[]){
    clock_t t1, t2, t3;
    t1 = clock();

    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }
    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//

    // ******************* Read Simulation Parameters *********//
    double_t re=1.0, s=7.0;
    double_t percentStrain = 15;
    double_t initialSearchRad = 1.0, finalSearchRad = 1.2;

    // Read input parameters from miscInp.dat file
    InputParameters miscInp = OPS::readKeyValueInput("miscInp.dat");
    re = std::stod( miscInp["re"] );
    percentStrain = std::stod( miscInp["percentStrain"] );
    initialSearchRad = std::stod( miscInp["initialSearchRad"] );
    finalSearchRad = std::stod( miscInp["finalSearchRad"] );

    s = (100 / (re*percentStrain))*log(2.0);

    std::ifstream coolFile("schedule.dat");
    assert(coolFile);
    std::vector<double_t> coolVec;
    double_t currGamma;

    while (coolFile >> currGamma) {
        coolVec.push_back(currGamma);
    }
    coolFile.close();

    // **********************************************************//

    // ***************** Create Bodies and Model ****************//
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

    // Prepare memory for energy, force
    double_t f;
    Eigen::VectorXd x(6*N), g(6*N);
    g.setZero(g.size());
    x.setZero(x.size());

    // Fill x with coords and rotVecs
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N), xrot(&(x(3*N)),3,N);
    xpos = coords;
    xrot = rotVecs;

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> posGrad(g.data(),3,N), rotGrad(&g(3*N),3,N);
    OPSMesh ops(N,f,xpos,xrot,posGrad,rotGrad);
    ops.setMorseDistance(re);
    s = 100*log(2.0)/(re*percentStrain);
    ops.setMorseWellWidth(s);
    ops.setSearchRadius(initialSearchRad);

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    // ****************************************************************//

    // ***************** Prepare Output Data file *********************//
    std::string fname = inputFileName.substr(0, inputFileName.find("."));
    std::stringstream sstm;

    std::string outputFileName;
    sstm << fname << "-Output.dat";
    outputFileName = sstm.str();
    sstm.str("");
    sstm.clear();

    ofstream outputFile;
    outputFile.open(outputFileName.c_str());
    outputFile << "#Step" << "\t"
               << "Gamma" << "\t"
               << "Asphericity" << "\t"
               << "Radius" << "\t"
               << "Volume" << "\t"
               << "Area" << "\t"
               << "MorseEnergy" << "\t"
               << "NormalityEn" << "\t"
               << "CircularityEn" << "\t"
               << "TotalOPSEnergy" << "\t"
               << "RMSAngleDeficit"
               << std::endl;
    // ******************************************************************//

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e7;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Calculate Average Edge Length
    double_t avgEdgeLen = ops.getAverageEdgeLength();
    std::cout << "Initial Avg Edge Length = " << avgEdgeLen << std::endl;

    // Renormalize positions such that avgEdgeLen = 1.0
    for(auto i=0; i < N; ++i){
        xpos.col(i) = xpos.col(i)/avgEdgeLen;
    }
    ops.setSearchRadius(finalSearchRad);
    ops.updatePolyData();
    ops.updateNeighbors();
    // ******************************************************************//

    t3 = clock();
    // ************************ SOLUTION LOOP **********************//
    int step=0;
    for(auto gamma : coolVec){
        std::cout<<"Iteration number = "<< step << std::endl;

        // Update OPS params
        ops.setFVK(gamma);

        // Solve the unconstrained minimization
        solver.solve();

        //Update polyData
        ops.updatePolyData();

        outputFile << step << "\t"
                   << gamma << "\t"
                   << ops.getAsphericity() << "\t"
                   << ops.getAverageRadius() << "\t"
                   << ops.getVolume() << "\t"
                   << ops.getArea() << "\t"
                   << ops.getMorseEnergy() << "\t"
                   << ops.getNormalityEnergy() << "\t"
                   << ops.getCircularityEnergy() << "\t"
                   << f << "\t"
                   << ops.getRMSAngleDeficit()
                   << std::endl;

        //********** Print relaxed configuration ************//
        sstm << fname << "-relaxed-" << step++ <<".vtk";
        std::string rName = sstm.str();
        ops.printVTKFile(rName);
        sstm.str("");
        sstm.clear();
    }
    // *****************************************************************************//
    t3 = clock() - t3;
    std::cout<<"Time for loop only = " << ((float)t3)/CLOCKS_PER_SEC << std::endl;

    outputFile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 1;
}
