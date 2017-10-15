#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>
#include "ALConstraint.h"
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
    std::string constraintType("NULL");
    enum Constraint{ AvgArea, AvgVol, ExactArea, ExactVol };
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
            >> temp >> searchRadFactor
            >> temp >> constraintType;

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
    else{
        std::cout<< "Invalid constraint type specified." << std::endl;
        exit(EXIT_FAILURE);
    }

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
    Eigen::Matrix3Xd coords(3,N);
    for(size_t i = 0; i < N; ++i){
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
    OPSParams params(D_e,re,s,b,initialSearchRad,gamma);
    OPSBody ops(N,f,xpos,xrot,posGrad,rotGrad,params);

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
                 << "Alpha" << "\t"
                 << "Beta" << "\t"
                 << "Gamma" << "\t"
                 << "Asphericity" << "\t"
                 << "Radius" << "\t"
                 << "Volume" << "\t"
                 << "MorseEnergy" << "\t"
                 << "NormalityEn" << "\t"
                 << "CircularityEn" << "\t"
                 << "TotalOPSEnergy" << "\t"
                 << "BrownianEnergy" << "\t"
                 << "ViscosityEnergy" << "\t"
                 << "TotalFunctional" <<"\t"
                 << "MSD" << "\t"
                 << "Neighbors" << "\t"
                 << "AvgEdgeLen" << "\t"
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

        // Update OPS params
        s = (100 / (avgEdgeLen*percentStrain))*log(2.0);
        params.updateParameter(OPSParams::gammaV, gamma);
        params.updateParameter(OPSParams::aV, s);

        // Solve the system without Brownian, Viscous and Volume constraints
        brown.setCoefficient(0.0);
        visco.setViscosity(0.0);
        constraint->setLagrangeCoeff(0.0);
        constraint->setPenaltyCoeff(0.0);
        std::cout<<"Solving the system at zero temperature..."<<std::endl;
        solver.solve();
        ops.updatePolyDataAndKdTree();
        ops.updateNeighbors();
        std::cout<<"Solving finished."<<std::endl;

        // Set up the constraint value as the zero temperature value
        double_t constrainedVal;
        if(type == AvgArea){
            double_t Ravg = ops.getAverageRadius();
             constrainedVal = 4*M_PI*Ravg*Ravg;
            std::cout<< "Constrained Area = " << constrainedVal << std::endl;
        }
        else if(type == AvgVol){
            double_t Ravg = ops.getAverageRadius();
            constrainedVal = 4*M_PI*Ravg*Ravg*Ravg/3;
            std::cout<< "Constrained Volume = " << constrainedVal << std::endl;
        }
        else if(type == ExactArea){
            constrainedVal = ops.getArea();
            std::cout<< "Constrained Area = " << constrainedVal << std::endl;
        }
        else if(type == ExactVol){
            constrainedVal = ops.getVolume();
            std::cout<< "Constrained Volume = " << constrainedVal << std::endl;
        }
        constraint->setConstraint(constrainedVal);

        // Set the viscosity and Brownian coefficient
        brownCoeff = beta*D_e/(alpha*avgEdgeLen);
        viscosity = brownCoeff/(alpha*avgEdgeLen);
        std::cout<< "Viscosity = " << viscosity << std::endl;
        std::cout<< "Brownian Coefficient = " << brownCoeff << std::endl;
        brown.setCoefficient(brownCoeff);
        visco.setViscosity(viscosity);

        // Update prevX
        prevX = x.head(3*N);

        //**************  INNER SOLUTION LOOP ******************//
        Eigen::Matrix3Xd averagePosition( 3, N ); /*!< Store avg particle positions */
        averagePosition = Eigen::Matrix3Xd::Zero(3,N);

        for (int viter = 0; viter < viterMax; viter++) {
            std::cout << std::endl
                      << "VISCOUS ITERATION: " << viter + stepCount
                      << std::endl
                      << std::endl;

            // Generate Brownian Kicks
            brown.generateParallelKicks();

            // Store data for Kabsch
            ops.updateDataForKabsch();

            // Set the starting guess for Lambda and K for Augmented Lagrangian
            constraint->setLagrangeCoeff(10.0);
            constraint->setPenaltyCoeff(1000.0);

            // *************** Augmented Lagrangian Loop ************** //
            double_t areaDiff = 1.0, areaTol = 1e-10;
            size_t alIter = 0, alMaxIter = 10;

            while( (areaDiff > areaTol) && (alIter < alMaxIter)){
                std::cout<< "Augmented Lagrangian iteration: " << alIter
                         << std::endl;

                // Solve the unconstrained minimization
                solver.solve();

                //Uzawa update
                constraint->uzawaUpdate();

                // Update termination check quantities
                alIter++;
                areaDiff = std::abs(constraint->getConstraintValue());

            }
            std::cout<< "areaDiff = " << areaDiff << "\talIter = "<< alIter
                     << std::endl << std::endl;
            // *********************************************************//

            // Apply Kabsch Algorithm
            ops.applyKabschAlgorithm();

            //Update kdTree, polyData and neighbors
            ops.updatePolyDataAndKdTree();
            ops.updateNeighbors();

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
                          << brown.getBrownianEnergy() << "\t"
                          << visco.getViscosityEnergy() << "\t"
                          << f << "\t"
                          << ops.getMeanSquaredDisplacement() << "\t"
                          << ops.getAverageNumberOfNeighbors() << "\t"
                          << avgEdgeLen
                          << std::endl;

            // Update prevX
            prevX = x.head(3*N);
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
    delete constraint;
    return 1;
}


