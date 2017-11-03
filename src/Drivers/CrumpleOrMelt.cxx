#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>
#include "ALConstraint.h"
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"

int main(int argc, char* argv[]){
    clock_t t1, t2, t3;
    t1 = clock();

    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }
    bool loggingOn = false;
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
    double_t De=1.0, re=1.0, s=7.0;
    double_t alpha=1.0, gamma=1.0;
    double_t percentStrain = 15, alpha_start=0.001, alpha_increment=0.001;
    double_t crumpledAsphericity = 0.01;

    std::string constraintType("NULL");
    enum Constraint{ AvgArea, AvgVol, ExactArea, ExactVol, ExactAreaAndVolume};
    //int lat_res=100, long_res=101;
    size_t viterMax = 500;
    size_t nameSuffix = 0;

    std::ifstream miscInpFile("miscInp.dat");
    assert(miscInpFile);
    std::string temp;
    miscInpFile
            >> temp >> De
            >> temp >> re
            >> temp >> constraintType
            >> temp >> alpha_start
            >> temp >> alpha_increment;

    miscInpFile.close();
    s = (100 / (re*percentStrain))*log(2.0);

    //Validate constraint type
    Constraint type;
    if(constraintType.compare("AverageArea") == 0){
        type = AvgArea;
    }
    else if(constraintType.compare("AverageVolume") == 0){
        type = AvgArea;
    }
    else if(constraintType.compare("ExactArea") == 0){
        type = ExactArea;
    }
    else if(constraintType.compare("ExactVolume") == 0){
        type = ExactVol;
    }
    else if(constraintType.compare("ExactAreaAndVolume") == 0){
        type = ExactAreaAndVolume;
    }
    else{
        std::cout<< "Invalid constraint type specified." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Input file should contain 3 columns
    // Gamma AreaConstrained PrintStep
    std::ifstream coolFile("cooling.dat");
    assert(coolFile);
    std::vector<std::vector<double_t> > coolVec;
    double_t currGamma, currPercentStrain, currArea, currPrintStep;

    std::string headerline;
    std::getline(coolFile, headerline);

    while (coolFile >> currGamma >> currPercentStrain >> currArea
           >> currPrintStep) {
        std::vector<double> currLine;
        currLine.push_back(currGamma);
        currLine.push_back(currPercentStrain);
        currLine.push_back(currArea);
        currLine.push_back(currPrintStep);
        coolVec.push_back(currLine);
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
    Eigen::VectorXd x(6*N), g(6*N), prevX(3*N);
    g.setZero(g.size());
    x.setZero(x.size());

    // Fill x with coords and rotVecs
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N), xrot(&(x(3*N)),3,N);
    xpos = coords;
    xrot = rotVecs;
    prevX = x.head(3*N);

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> posGrad(g.data(),3,N), rotGrad(&g(3*N),3,N);
    OPSMesh ops(N,f,xpos,xrot,posGrad,rotGrad);
    ops.setMorseDistance(re);
    ops.setMorseEnergy(De);
    s = 100*log(2.0)/(re*percentStrain);
    ops.setMorseWellWidth(s);

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    double_t brownCoeff = 1.0, viscosity = 1.0;
    BrownianBody brown(3*N,brownCoeff,f,thermalX,thermalG,prevX);
    ViscosityBody visco(3*N,viscosity,f,thermalX,thermalG,prevX);

    // Create the Augmented Lagrangian volume constraint body
    ALConstraint* constraint;
    if(type == AvgArea){
        constraint = new AvgAreaConstraint(N, f, xpos, posGrad);
    }
    else if(type == AvgVol){
        constraint = new AvgVolConstraint(N, f, xpos, posGrad);
    }
    else if(type == ExactArea){
        vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
        constraint = new ExactAreaConstraint(N, f, xpos, posGrad, poly);
    }
    else if(type == ExactVol){
        vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
        constraint = new ExactVolConstraint(N, f, xpos, posGrad, poly);
    }
    else if(type == ExactAreaAndVolume){
        vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
        constraint = new ExactAreaVolConstraint(N, f, xpos, posGrad, poly);
        constraint->setTolerance(1e-8);
    }

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);
    model.addBody(constraint);
    // ****************************************************************//

    // ***************** Prepare Output Data file *********************//
    std::string fname = inputFileName.substr(0, inputFileName.find("."));
    std::stringstream sstm;

    std::string dataOutputFile;
    sstm << fname << "-DetailOutput.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();

    ofstream innerLoopFile;
    innerLoopFile.open(dataOutputFile.c_str());
    innerLoopFile << "#Step" << "\t"
                  <<"ParaviewStep" << "\t"
                 << "Alpha" <<"\t"
                 << "Gamma" << "\t"
                 << "Asphericity" << "\t"
                 << "Radius" << "\t"
                 << "Volume" << "\t"
                 << "Area" << "\t"
                 << "MorseEnergy" << "\t"
                 << "NormalityEn" << "\t"
                 << "CircularityEn" << "\t"
                 << "TotalOPSEnergy" << "\t"
                 << "BrownianEnergy" << "\t"
                 << "ViscosityEnergy" << "\t"
                 << "TotalFunctional" <<"\t"
                 << "MSD"
                 << std::endl;

    ofstream outerLoopFile;
    sstm << fname << "-CrumplingOutput.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    outerLoopFile.open(dataOutputFile.c_str());
    outerLoopFile << "#BigStep" <<"\t"
                  << "Gamma" << "\t"
                  << "PercentStrain" << "\t"
                  << "CrumplingAlpha" <<"\t"
                  << "CrumplingRadius"  <<"\t"
                  << "Asphericity" << "\t"
                  << "Volume"
                  << std::endl;
    // ******************************************************************//

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e7;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    solver.turnOffLogging();
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Calculate Average Edge Length
    double_t avgEdgeLen = ops.getAverageEdgeLength();
    if(loggingOn)
        std::cout << "Initial Avg Edge Length = " << avgEdgeLen << std::endl;

    // Renormalize positions such that avgEdgeLen = 1.0
    for(auto i=0; i < N; ++i){
        xpos.col(i) = xpos.col(i)/avgEdgeLen;
    }

    // Update the OPSBody member variables as per new positions
    ops.updatePolyData();
    ops.updateNeighbors();
    ops.saveInitialPosition(); /*!< For Mean Squared Displacement or reset*/
    avgEdgeLen = ops.getAverageEdgeLength();
    if(loggingOn)
        std::cout << "After renormalizing, Avg Edge Length = "
              << avgEdgeLen << std::endl;
    // ******************************************************************//

    t3 = clock();
    // ************************ OUTER SOLUTION LOOP **********************//
    int printStep, step=0, z = 0;
    int paraviewStep = -1;
    for(auto line : coolVec){
        gamma = line[0];
        percentStrain = line[1];
        double_t constrainedVal = line[2];
        printStep = (int)line[3];

        // Update OPS params
        s = (100 / (avgEdgeLen*percentStrain))*log(2.0);
        ops.setFVK(gamma);
        ops.setMorseWellWidth(s);

        // Set up the area constraint value
        constraint->setConstraint(constrainedVal);

        // Update prevX
        prevX = x.head(3*N);

        // ************************* CRUMPLING LOOP ***************************//
        bool crumpled = false;
        alpha = alpha_start - alpha_increment;
        double_t asphericity = 0.0;

        while(!crumpled){

            // Increment alpha
            alpha += alpha_increment;

            // Set the viscosity and Brownian coefficient
            brownCoeff = De/(alpha*avgEdgeLen);
            viscosity = brownCoeff/(alpha*avgEdgeLen);
            if(loggingOn){
                std::cout<< "Viscosity = " << viscosity << std::endl;
                std::cout<< "Brownian Coefficient = " << brownCoeff << std::endl;
            }
            brown.setCoefficient(brownCoeff);
            visco.setViscosity(viscosity);

            //**************  INNER LOOP ******************//
            for (auto viter = 0; viter < viterMax; viter++) {

                // Generate Brownian Kicks
                brown.generateParallelKicks();

                // Store data for Kabsch
                ops.updateDataForKabsch();

                // Set the starting guess for Lambda and K for constraints
                constraint->setLagrangeCoeff(10.0);
                constraint->setPenaltyCoeff(1000.0);

                // *************** Augmented Lagrangian Loop ************** //
                bool constraintMet = false;
                size_t alIter = 0, alMaxIter = 10;

                while( !constraintMet && (alIter < alMaxIter)){
                    if(loggingOn)
                        std::cout<< "Augmented Lagrangian iteration: " << alIter
                             << std::endl;

                    // Solve the unconstrained minimization
                    solver.solve();

                    //Uzawa update
                    constraint->uzawaUpdate();

                    // Update termination check quantities
                    alIter++;
                    constraintMet = constraint->constraintSatisfied();
                }
                if(loggingOn){
                    constraint->printCompletion();
                    std::cout<< "Constraint satisfied in "<< alIter << " iterations."
                         << std::endl << std::endl;
                }
                // *********************************************************//

                // Apply Kabsch Algorithm
                ops.applyKabschAlgorithm();

                //Update kdTree, polyData and neighbors
                ops.updatePolyData();
                ops.updateNeighbors();

                //********** Print relaxed configuration ************//
                //We will print only after every currPrintStep iterations
                if (viter % printStep == 0) {
                    paraviewStep++;
                    sstm << fname << "-relaxed-" << nameSuffix++ <<".vtk";
                    std::string rName = sstm.str();
                    ops.printVTKFile(rName);
                    sstm.str("");
                    sstm.clear();
                }

                int paraviewStepPrint;
                paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;
                asphericity = ops.getAsphericity();

                innerLoopFile << step++ << "\t"
                              << paraviewStepPrint <<"\t"
                              << alpha <<"\t"
                              << gamma << "\t"
                              << asphericity << "\t"
                              << ops.getAverageRadius() << "\t"
                              << ops.getVolume() << "\t"
                              << ops.getArea() << "\t"
                              << ops.getMorseEnergy() << "\t"
                              << ops.getNormalityEnergy() << "\t"
                              << ops.getCircularityEnergy() << "\t"
                              << ops.getTotalEnergy() << "\t"
                              << brown.getBrownianEnergy() << "\t"
                              << visco.getViscosityEnergy() << "\t"
                              << f << "\t"
                              << ops.getMeanSquaredDisplacement()
                              << std::endl;

                if( asphericity >= crumpledAsphericity ){
                    crumpled = true;
                    break;
                }

                // Update prevX
                prevX = x.head(3*N);
            }
            // ********************************************************//
        }
        //******************************************************************//

        // Write the crumpling alpha values and statistics to output file
        outerLoopFile << z++ <<"\t"
                      << gamma << "\t"
                      << percentStrain << "\t"
                      << alpha << "\t"
                      << ops.getAverageRadius() <<"\t"
                      << asphericity <<"\t"
                      << ops.getVolume()
                      << std::endl;

        // Reset OPSMesh for next gamma value
        ops.resetToInitialPositions();
    }
    // *****************************************************************************//

    innerLoopFile.close();
    outerLoopFile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    delete constraint;
    return 1;
}


