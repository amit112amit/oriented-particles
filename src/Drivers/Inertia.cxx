#include <cstdlib>
#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>
#include <Eigen/Eigenvalues>
#include "ALConstraint.h"
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"

using namespace OPS;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Vector3d Vector3d;

int main(int argc, char* argv[]){
    clock_t t1, t2, t3;
    t1 = clock();

    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }

    //******************** Optional parameters ********************//
    bool loggingOn = false;

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];
    std::stringstream taskId;
    std::stringstream sstm;
    taskId << std::getenv("SGE_TASK_ID");
    std::cout << "SGE_TASK_ID = " << taskId.str() << std::endl;

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->ReadAllVectorsOn();
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//

    // ******************* Read Simulation Parameters *********//

    // This flag determines if its a new simulation or a continuation

    double_t re=1.0, s=7.0;
    double_t alpha=1.0, beta=1.0, gamma=1.0;
    double_t percentStrain = 15;

    enum Constraint{ AvgArea, AvgVol, ExactArea, ExactVol, ExactAreaAndVolume};
    std::string constraintType("NULL"), baseFileName;
    size_t viterMax = 1000;
    size_t nameSuffix = 0;
    size_t step = 0;

    InputParameters miscInp = OPS::readKeyValueInput( "miscInp.dat" );
    re = std::stod( miscInp["re"] );
    constraintType = miscInp["constraintType"];
    baseFileName = miscInp["baseFileName"];

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

    // Input file should contain the following columns
    // Alpha Beta Gamma PercentStrain AreaConstraint NumIterations PrintStep
    std::ifstream coolFile("schedule.dat");
    assert(coolFile);
    std::vector<std::vector<double_t> > coolVec;
    double_t currAlpha, currBeta, currGamma, currPercentStrain,
                    currViterMax, currPrintStep, currArea;

    std::string headerline;
    std::getline(coolFile, headerline);

    while (coolFile >> currAlpha >> currBeta >> currGamma >> currPercentStrain
                    >> currArea >> currViterMax >> currPrintStep) {
        std::vector<double> currLine;
        currLine.push_back(currAlpha);
        currLine.push_back(currBeta);
        currLine.push_back(currGamma);
        currLine.push_back(currPercentStrain);
        currLine.push_back(currArea);
        currLine.push_back(currViterMax);
        currLine.push_back(currPrintStep);
        coolVec.push_back(currLine);
    }
    coolFile.close();

    // **********************************************************//

    // ***************** Create Bodies and Model ****************//
    // Set number of OPS particles
    size_t N = mesh->GetNumberOfPoints();

    // Read point coordinates from input mesh
    Eigen::Matrix3Xd coords(3,N);
    for(auto i = 0; i < N; ++i){
        Eigen::Vector3d cp = Eigen::Vector3d::Zero();
        mesh->GetPoint(i, &(cp(0)));
        coords.col(i) = cp;
    }

    // Generate initial rotation vectors either from starting point coordinates
    Eigen::Matrix3Xd rotVecs(3,N);
    // Generate rotation vectors from input point coordinates
    OPSBody::initialRotationVector(coords, rotVecs);

    // Prepare memory for energy, force
    double_t f;
    Eigen::VectorXd x(6*N), g(6*N), prevX(3*N);
    g.setZero(g.size());
    x.setZero(x.size());

    // Fill x with coords and rotVecs
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N), xrot(&(x(3*N)),3,N),
            prevPos(prevX.data(),3,N);
    xpos = coords;
    xrot = rotVecs;
    prevX = x.head(3*N);

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> posGrad(g.data(),3,N), rotGrad(&g(3*N),3,N);
    OPSMesh ops(N,f,xpos,xrot,posGrad,rotGrad,prevPos);
    ops.setMorseDistance(re);
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

    // ***************** Prepare Output Data files *********************//
    // Identify the Input structure name
    std::string fname = baseFileName;
    std::string dataOutputFile;
    sstm << fname << "-DetailedOutput-" << taskId.str() << ".dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    // Detailed output data file
    ofstream detailedOP;
    detailedOP.open(dataOutputFile.c_str(), std::ofstream::out);
    detailedOP
                    << "Gamma" << "\t"
                    << "Beta" << "\t"
                    << "Asphericity" << "\t"
                    << "MorseEn"  <<"\t"
                    << "NormEn"  <<"\t"
                    << "CircEn"  <<"\t"
                    << "BrownEn"  <<"\t"
                    << "ViscoEn"  <<"\t"
                    << "MSD" << "\t"
                    << "RMSAngleDeficit" << "\t"
                    << "I1" << "\t"
                    << "I2" << "\t"
                    << "I3"
                    << std::endl;

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e5;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    solver.turnOffLogging();
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Calculate Average Edge Length
    double_t avgEdgeLen = ops.getAverageEdgeLength();
    if(loggingOn){
        std::cout << "Initial Avg Edge Length = " << avgEdgeLen
                  << std::endl;
    }
    // Renormalize positions such that avgEdgeLen = 1.0
    for(auto i=0; i < N; ++i){
        xpos.col(i) = xpos.col(i)/avgEdgeLen;
    }

    // Update the OPSBody member variables as per new positions
    ops.updatePolyData();
    ops.updateNeighbors();
    ops.saveInitialPosition(); /*!< For Mean Squared Displacement */
    avgEdgeLen = ops.getAverageEdgeLength();
    if(loggingOn)
        std::cout << "After renormalizing, Avg Edge Length = "
                  << avgEdgeLen << std::endl;
    // ******************************************************************//

    //Create an eigenvalue solver for inertia tensor
    Eigen::SelfAdjointEigenSolver< Matrix3d > saes;

    // ************************ OUTER SOLUTION LOOP **********************//
    size_t printStep;
    for(int z=0; z < coolVec.size(); z++){
        alpha = coolVec[z][0];
        beta = coolVec[z][1];
        gamma = coolVec[z][2];
        percentStrain = coolVec[z][3];
        double_t constrainedVal = coolVec[z][4];
        viterMax = coolVec[z][5];
        printStep = (int)coolVec[z][6];

        // Update OPS params
        s = (100 / percentStrain)*log(2.0);
        ops.setFVK(gamma);
        ops.setMorseWellWidth(s);

        // Set up the constraint value as the zero temperature value
        constraint->setConstraint(constrainedVal);

        // For the very first iteration solve at zero temperature first
        if( z == 0 ){
            brown.setCoefficient(0.0);
            visco.setViscosity(0.0);
            solver.solve();
        }

        // Update prevX
        prevX = x.head(3*N);

        // Set the viscosity and Brownian coefficient
        viscosity = alpha;
        brownCoeff = std::sqrt( 2*alpha/beta );
        if(loggingOn){
            std::cout<< "Viscosity = " << viscosity << std::endl;
            std::cout<< "Brownian Coefficient = " << brownCoeff
                     << std::endl;
        }
        brown.setCoefficient(brownCoeff);
        visco.setViscosity(viscosity);

        //**************  INNER SOLUTION LOOP ******************//
        // Average energy across time steps
        double_t avgTotalEnergy = 0.0;

        for (int viter = 0; viter < viterMax; viter++) {
            if(loggingOn)
                std::cout << std::endl
                          << "VISCOUS ITERATION: " << step
                          << std::endl
                          << std::endl;

            // Generate Brownian Kicks
            brown.generateParallelKicks();

            // Set the starting guess for Lambda and K for
            // Augmented Lagrangian
            constraint->setLagrangeCoeff(10.0);
            constraint->setPenaltyCoeff(1000.0);

            // *************** Augmented Lagrangian Loop ************** //
            bool constraintMet = false;
            size_t alIter = 0, alMaxIter = 10;

            while( !constraintMet && (alIter < alMaxIter)){
                if(loggingOn)
                    std::cout<<"Augmented Lagrangian iteration: "
                            << alIter
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
                std::cout<< "Constraint satisfied in "<< alIter
                         << " iterations."
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
            if (viter % printStep == 0 && printStep <= viterMax) {
                sstm << fname << "-" << taskId.str() << "-relaxed-"
                     << nameSuffix++ <<".vtk";
                std::string rName = sstm.str();
                ops.printVTKFile(rName);
                sstm.str("");
                sstm.clear();
            }

            // Calculate the inertia tensor and its eigenvalues
            Matrix3d M = Matrix3d::Zero();
            for( auto ptId = 0; ptId < N; ++ptId ){
                Vector3d xk = xpos.col(ptId);
                M += xk.dot(xk)*Matrix3d::Identity() - xk*xk.transpose();
            }
            // Now calculate the eigen values
            saes.compute( M, Eigen::EigenvaluesOnly );
            Vector3d eigV = saes.eigenvalues();

            double_t I1 = eigV[0];
            double_t I2 = eigV[1];
            double_t I3 = eigV[2];

            std::vector<double_t> msds(2,0);
            msds = ops.getMSD();

            // Write output to data file
            detailedOP
                            << gamma << "\t"
                            << beta << "\t"
                            << ops.getAsphericity() << "\t"
                            << ops.getMorseEnergy() << "\t"
                            << ops.getNormalityEnergy() << "\t"
                            << ops.getCircularityEnergy() << "\t"
                            << brown.getBrownianEnergy() << "\t"
                            << visco.getViscosityEnergy() << "\t"
                            << msds[0] << "\t"
                            << ops.getRMSAngleDeficit() << "\t"
                            << I1 << "\t"
                            << I2 << "\t"
                            << I3
                            << std::endl;

            // Update prevX
            prevX = x.head(3*N);

            if( loggingOn ){
                avgTotalEnergy=(avgTotalEnergy*step+f)/(step+1);
                std::cout<< " Average Total Energy = "
                         << avgTotalEnergy << std::endl;
            }
            step++;
        }
        //************************************************//

    }
    // **********************************************************************//
    detailedOP.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    delete constraint;
    return 1;
}
