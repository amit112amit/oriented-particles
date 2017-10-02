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
    double_t D_e=1.0, re=1.0, s=7.0, b=1.0;
    double_t alpha=1.0, beta=1.0, gamma=1.0;
    double_t percentStrain = 15;
    double_t initialSearchRad = 1.0, finalSearchRad = 1.2,
            searchRadFactor = 1.3;
    //int lat_res=100, long_res=101;
    size_t viterMax = 1000;
    size_t nameSuffix = 0;

    std::ifstream miscInpFile("miscInp.dat");
    assert(miscInpFile);
    std::string temp;
    miscInpFile
            >> temp >> D_e
            >> temp >> re
            >> temp >> b
            >> temp >> initialSearchRad
            >> temp >> finalSearchRad
            >> temp >> searchRadFactor;
    miscInpFile.close();
    s = (100 / (re*percentStrain))*log(2.0);

    std::ifstream coolFile("cooling.dat");
    assert(coolFile);
    std::vector<std::vector<double_t> > coolVec;
    double_t currAlpha, currBeta, currGamma, currPercentStrain,
            currViterMax, currPrintStep;

    std::string headerline;
    std::getline(coolFile, headerline);

    while (coolFile >> currAlpha >> currBeta >> currGamma >>
           currPercentStrain >> currViterMax >> currPrintStep) {
        std::vector<double> currLine;
        currLine.push_back(currAlpha);
        currLine.push_back(currBeta);
        currLine.push_back(currGamma);
        currLine.push_back(currPercentStrain);
        currLine.push_back(currViterMax);
        currLine.push_back(currPrintStep);
        coolVec.push_back(currLine);
    }
    coolFile.close();

    // **********************************************************//

    // ***************** Create Bodies and Model ****************//
    // Set number of OPS particles
    size_t N = mesh->GetNumberOfPoints();

    // Generate Rotation Vectors from input point coordinates
    Eigen::Map<Eigen::Matrix3Xd> coords(
                (double_t*)mesh->GetPoints()->GetData()->GetVoidPointer(0),
                3,N);
    Eigen::Matrix3Xd rotVecs(3,N);
    OPSBody::initialRotationVector(coords, rotVecs);

    // Prepare memory for energy and force
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
    OPSParams params(D_e,re,s,b,initialSearchRad,gamma);
    OPSBody ops(N,f,xpos,xrot,posGrad,rotGrad,params);

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    double brownCoeff = 1.0, viscosity = 1.0;
    BrownianBody brown(3*N,brownCoeff,f,thermalX,thermalG);
    ViscosityBody visco(3*N,viscosity,f,thermalX,thermalG);

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);
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
                 << "Alpha" << "\t"
                 << "Beta" << "\t"
                 << "Gamma" << "\t"
                 << "Asphericity" << "\t"
                 << "Radius" << "\t"
                 << "Volume" << "\t"
                 << "MorseEnergy" << "\t"
                 << "NormalityEn" << "\t"
                 << "CircularityEn" << "\t"
                 << "BrownianEnergy" << "\t"
                 << "ViscosityEnergy" << "\t"
                 << "PressureVolEn" << "\t"
                 << "TotalFunctional" <<"\t"
                 << "MSD" << "\t"
                 << "Neighbors"
                 << std::endl;

    ofstream outerLoopFile;
    sstm << fname << "-AverageOutput.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    outerLoopFile.open(dataOutputFile.c_str());
    outerLoopFile << "#BigStep" <<"\t"
                  << "PercentStrain" <<"\t"
                  << "Alpha" << "\t"
                  << "Beta" << "\t"
                  << "Gamma" << "\t"
                  << "PercentStrain" << "\t"
                  << "Radius"  <<"\t"
                  << "Asphericity"
                  << std::endl;
    // ******************************************************************//

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 10000;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Calculate Average Edge Length
    double_t avgEdgeLen = ops.getAverageEdgeLength();
    std::cout << "Initial Avg Edge Length = " << avgEdgeLen << std::endl;

    // Renormalize positions such that avgEdgeLen = 1.0
    for(size_t i=0; i < N; ++i){
        xpos.col(i) = xpos.col(i)/avgEdgeLen;
    }
    params.updateParameter(OPSParams::searchRadiusV, finalSearchRad);
    avgEdgeLen = ops.getAverageEdgeLength();
    std::cout << "After renormalizing, Avg Edge Length = "
              << avgEdgeLen << std::endl;

    // Update the OPSBody member variables as per new positions
    params.updateParameter(OPSParams::searchRadiusV,
                           searchRadFactor*avgEdgeLen);
    ops.updatePolyDataAndKdTree();
    ops.updateNeighbors();
    ops.saveInitialPosition(); /*!< For Mean Squared Displacement */
    // ******************************************************************//

    // ************************ OUTER SOLUTION LOOP **********************//
    int printStep, stepCount = 0, step=0;
    int paraviewStep = -1;
    for(int z=0; z < coolVec.size(); z++){
        alpha = coolVec[z][0];
        beta = coolVec[z][1];
        gamma = coolVec[z][2];
        percentStrain = coolVec[z][3];
        viterMax = coolVec[z][4];
        printStep = (int)coolVec[z][5];

        s = (100 / (avgEdgeLen*percentStrain))*log(2.0);
        brownCoeff = beta*D_e/(alpha*avgEdgeLen);
        viscosity = brownCoeff/(alpha*avgEdgeLen);
        std::cout<< "Viscosity = " << viscosity << std::endl;
        std::cout<< "Brownian Coefficient = " << brownCoeff << std::endl;
        brown.setCoefficient(brownCoeff);
        visco.setViscosity(viscosity);

        params.updateParameter(OPSParams::gammaV, gamma);
        params.updateParameter(OPSParams::aV, s);

        //**************  INNER SOLUTION LOOP ******************//
        Eigen::Matrix3Xd averagePosition( 3, N ); /*!< Store avg particle positions */
        averagePosition = Eigen::Matrix3Xd::Zero(3,N);

        for (int viter = 0; viter < viterMax; viter++) {
            //Ensure only next nearest neighbor interactions
            avgEdgeLen = ops.getAverageEdgeLength();
            params.updateParameter(OPSParams::searchRadiusV,
                                   searchRadFactor*avgEdgeLen);
            std::cout << std::endl
                 << "VISCOUS ITERATION: " << viter + stepCount
                 << std::endl
                 << std::endl;

            // Store the curren position for applying Kabsch algorithm
            ops.updateDataForKabsch();

            // Generate Brownian Kicks
            brown.generateParallelKicks();

            // Solve the model
            solver.solve();

            // Apply Kabsch Algorithm
            ops.applyKabschAlgorithm();

            // Add current solution to average position data
            averagePosition += xpos;

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

            innerLoopFile << step++ << "\t"
                          << paraviewStepPrint <<"\t"
                          << alpha <<"\t"
                          << beta << "\t"
                          << gamma << "\t"
                          << ops.getAsphericity() << "\t"
                          << ops.getAverageRadius() << "\t"
                          << ops.getVolume() << "\t"
                          << ops.getMorseEnergy() << "\t"
                          << ops.getNormalityEnergy() << "\t"
                          << ops.getCircularityEnergy() << "\t"
                          << ops.getTotalEnergy() << "\t"
                          << visco.getViscosityEnergy() << "\t"
                          << f << "\t"
                          << ops.getMeanSquaredDisplacement() << "\t"
                          << ops.getAverageNumberOfNeighbors()
                          << std::endl;

            // Update viscosity body reference
            visco.viscousStep();
        }
        //************************************************//

        // Calculate the average particle positions and avg radius
        double avgShapeRad = 0.0, avgShapeAsph = 0.0;
        vtkSmartPointer<vtkPoints> avgPos =
                vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> avgPosData =
                vtkSmartPointer<vtkDoubleArray>::New();
        averagePosition = averagePosition / viterMax;
        avgShapeRad = averagePosition.colwise().norm().sum()/N;
        void *avgPosPtr = (void*)averagePosition.data();
        avgPosData->SetVoidArray(avgPosPtr,3*N,1);
        avgPosData->SetNumberOfComponents(3);
        avgPos->SetData(avgPosData);

        sstm << fname << "-AvgShape-"<< z << ".vtk";
        dataOutputFile = sstm.str();
        sstm.str("");
        sstm.clear();

        // Print the average shape
        delaunay3DSurf(avgPos, dataOutputFile);

        // Calculate the asphericity of the average shape        
        Eigen::RowVectorXd R(N);
        R = averagePosition.colwise().norm();
        avgShapeAsph += ((R.array() - avgShapeRad).square()).sum();
        avgShapeAsph /= (N*avgShapeRad*avgShapeRad);

        outerLoopFile << z <<"\t"
                      << percentStrain << "\t"
                      << alpha << "\t"
                      << beta << "\t"
                      << gamma << "\t"
                      << percentStrain << "\t"
                      << avgShapeRad  << "\t"
                      << avgShapeAsph
                      << std::endl;

    }
    // *****************************************************************************//

    innerLoopFile.close();
    outerLoopFile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;

    return 1;
}
