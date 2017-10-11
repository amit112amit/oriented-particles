#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSBody.h"
#include "ViscosityBody.h"

int main(int argc, char* argv[]){

    clock_t t1, t2;
    t1 = clock();

    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }
    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];

    vtkSmartPointer<vtkPolyDataReader> reader =
            vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//

    // ******************* Read Simulation Parameters *********//
    double_t alpha=1.0, beta=1.0;
    size_t viterMax = 1000;
    size_t nameSuffix = 0;

    std::ifstream coolFile("cooling.dat");
    assert(coolFile);
    std::vector<std::vector<double_t> > coolVec;
    double_t currAlpha, currBeta, currViterMax, currPrintStep;

    std::string headerline;
    std::getline(coolFile, headerline);

    while (coolFile >> currAlpha >> currBeta >> currViterMax
           >> currPrintStep) {
        std::vector<double> currLine;
        currLine.push_back(currAlpha);
        currLine.push_back(currBeta);
        currLine.push_back(currViterMax);
        currLine.push_back(currPrintStep);
        coolVec.push_back(currLine);
    }
    coolFile.close();

    // **********************************************************//

    // ***************** Create Bodies and Model ****************//
    // Set number of OPS particles
    size_t N = mesh->GetNumberOfPoints();

    // Prepare memory for energy and force
    double_t f;
    Eigen::VectorXd x(3*N), g(3*N), prevX(3*N);
    g.setZero(g.size());
    x.setZero(x.size());
    prevX.setZero(prevX.size());

    // Fill x and prevX with coordinates
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N);
    for(size_t i = 0; i < N; ++i){
        Eigen::Vector3d cp = Eigen::Vector3d::Zero();
        mesh->GetPoint(i, &(cp(0)));
        xpos.col(i) = cp;
    }
    prevX = x;

    // Create Brownian and Viscosity bodies    
    double brownCoeff = 1.0, viscosity = 1.0;
    BrownianBody brown(3*N,brownCoeff,f,x,g,prevX);
    ViscosityBody visco(3*N,viscosity,f,x,g,prevX);

    // Create Model
    Model model(3*N,f,g);
    model.addBody(&brown);
    model.addBody(&visco);
    // ****************************************************************//

    // ***************** Prepare Output Data file *********************//
    std::string fname = inputFileName.substr(0, inputFileName.find("."));
    std::stringstream sstm;

    std::string dataOutputFile;
    sstm << fname << "-Output.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();

    ofstream innerLoopFile;
    innerLoopFile.open(dataOutputFile.c_str());
    innerLoopFile << "#Step" << "\t"
                  <<"ParaviewStep" << "\t"
                 << "Alpha" << "\t"
                 << "Beta" << "\t"
                 << "BrownianEnergy" << "\t"
                 << "ViscosityEnergy" << "\t"
                 << "TotalFunctional" <<"\t"
                 << std::endl;
    // ******************************************************************//

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 10000;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    // *****************************************************************//

    // ************************ SOLUTION LOOP **********************//
    int printStep, stepCount = 0, step=0;
    int paraviewStep = -1;
    for(int z=0; z < coolVec.size(); z++){
        alpha = coolVec[z][0];
        beta = coolVec[z][1];
        viterMax = coolVec[z][2];
        printStep = (int)coolVec[z][3];

        brownCoeff = beta/alpha;
        viscosity = brownCoeff/alpha;
        std::cout<< "Viscosity = " << viscosity << std::endl;
        std::cout<< "Brownian Coefficient = " << brownCoeff << std::endl;
        brown.setCoefficient(brownCoeff);
        visco.setViscosity(viscosity);

        //**************  INNER SOLUTION LOOP ******************//
        for (int viter = 0; viter < viterMax; viter++) {
            std::cout << std::endl
                 << "VISCOUS ITERATION: " << viter + stepCount
                 << std::endl
                 << std::endl;

            // Generate Brownian Kicks
            brown.generateParallelKicks();

            // Solve the model
            solver.solve();

            //********** Print relaxed configuration ************//
            //We will print only after every currPrintStep iterations
            if (viter % printStep == 0) {
                paraviewStep++;
                sstm << fname << "-relaxed-" << nameSuffix++ <<".vtk";
                std::string rName = sstm.str();
                brown.printVTKFile(rName);
                sstm.str("");
                sstm.clear();
            }

            int paraviewStepPrint;
            paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;

            innerLoopFile << step++ << "\t"
                          << paraviewStepPrint <<"\t"
                          << alpha <<"\t"
                          << beta << "\t"
                          << brown.getBrownianEnergy() << "\t"
                          << visco.getViscosityEnergy() << "\t"
                          << f << "\t"
                          << std::endl;

            // Store the current state for next iteration
            prevX = x;
        }
        //************************************************//

    }
    // *****************************************************************************//

    innerLoopFile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;

    return 1;
}

