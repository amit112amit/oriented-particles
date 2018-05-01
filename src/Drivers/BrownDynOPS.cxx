#include <cstdlib>
#include <random>
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
#include "HelperFunctions.h"

using namespace OPS;
using namespace std;

int main(int argc, char* argv[]){
    clock_t t1;
    t1 = clock();

    if (argc != 2) {
        cout << "usage: " << argv[0] << " <filename>\n";
        return -1;
    }

    //******************** Create variables ********************//
    std::string inFile = argv[1];
    std::string baseFileName = inFile.substr(0,inFile.length() - 4);
    std::stringstream taskId;
    std::stringstream sstm;

    double_t alpha=1.0, beta=1.0, gamma=1.0, re=1.0, f = 0,
                    percentStrain = 15,
                    s=(100/(re*percentStrain))*log(2.0);

    size_t viterMax, nameSuffix=0, step=0, N=10, saveFreq = 100;

    Eigen::VectorXd x(6*N), prevX(3*N), g(6*N);
    Eigen::Matrix3Xd coords(3,N), rotVecs(3,N), initPos(3,N);
    std::vector<vtkIdType> neighbors(N);
    std::mt19937 engine;
    auto rng = std::normal_distribution<double_t>(0.0,1.0);
    //**********************************************************//

    taskId << std::getenv("SGE_TASK_ID");
    std::cout << "SGE_TASK_ID = " << taskId.str() << std::endl;

    //******************** Resume or fresh start? ********************//
    // Create a SimulationState
    SimulationState state;

    // Check if a state file exists in current directory
    sstm << "SimulationState-" << taskId.str() << ".dat";
    string stateFileName = sstm.str();
    sstm.str("");
    sstm.clear();
    ifstream stateFile(stateFileName);

    if(stateFile.good()){
        // Read previously stored state to resume the simulation
        state = SimulationState::readFromFile(stateFileName);
        nameSuffix = state.getNameSuffix();
        nameSuffix = nameSuffix > 0? nameSuffix + 1: 0;
        step = state.getStep() + 1;
        N = state.getN();
        engine = state.getRandomEngine();
        rng = state.getRandomGenerator();
        // Resize matrices and vectors
        x.resize(6*N);
        g.resize(6*N);
        prevX.resize(3*N);
        initPos.resize(3,N);
        coords.resize(3,N);
        rotVecs.resize(3,N);
        neighbors.resize(N);
        // Read the vectors into local variables
        x = state.getX();
        prevX = state.getPrevX();
        neighbors = state.getNeighbors();
        initPos = state.getInitPos();
        // Copy position and rotation vectors into coords, rotVecs
        for( auto i=0; i < N; ++i ){
            size_t si = 3*i;
            coords.col(i) << x(si), x(si+1), x(si+2);
        }
        for( auto i=0; i < N; ++i ){
            size_t si = 3*(N+i);
            rotVecs.col(i) << x(si), x(si+1), x(si+2);
        }
    }
    else{
        // No state file found. So we need to start a new simulation
        vtkSmartPointer<vtkPolyData> mesh;
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(inFile.c_str());
        reader->ReadAllVectorsOn();
        reader->Update();
        mesh = reader->GetOutput();
        N = mesh->GetNumberOfPoints();
        // Resize matrices and vectors
        x.resize(6*N);
        g.resize(6*N);
        prevX.resize(3*N);
        initPos.resize(3,N);
        coords.resize(3,N);
        rotVecs.resize(3,N);
        neighbors.resize(N);
        // Read point coordinates from input mesh
        for(auto i = 0; i < N; ++i){
            Eigen::Vector3d cp = Eigen::Vector3d::Zero();
            mesh->GetPoint(i, &(cp(0)));
            coords.col(i) = cp;
        }
        // Renormalize by the average edge length
        double_t avgEdgeLen = getPointCloudAvgEdgeLen(inFile);
        for(auto i=0; i < N; ++i){
            coords.col(i) = coords.col(i)/avgEdgeLen;
        }
        // Generate rotation vectors from input point coordinates
        OPSBody::initialRotationVector(coords, rotVecs);
    }
    // ****************************************************************//

    // ***************** Create Bodies and Model ****************//

    // Fill x with coords and rotVecs
    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(),3,N), xrot(&(x(3*N)),3,N),
                    prevPos(prevX.data(),3,N);
    xpos = coords;
    xrot = rotVecs;

    // Create OPSBody
    g.setZero(g.size());
    Eigen::Map<Eigen::Matrix3Xd> posGrad(g.data(),3,N), rotGrad(&g(3*N),3,N);
    OPSMesh ops(N,f,xpos,xrot,posGrad,rotGrad,prevPos);
    ops.setMorseDistance(re);
    ops.setMorseWellWidth(s);
    ops.updatePolyData();
    ops.updateNeighbors();

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),3*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),3*N,1);
    double_t brownCoeff = 1.0, viscosity = 1.0;
    BrownianBody brown(3*N,brownCoeff,f,thermalX,thermalG,prevX);

    // Set ops history variables
    if( stateFile.good() ){
        ops.setInitialNeighbors(neighbors);
        ops.setInitialPositions(initPos);
        brown.setRandomEngine(engine);
        brown.setRandomGenerator(rng);
    }
    else{
        prevX = x.head(3*N);
        neighbors = ops.getInitialNeighbors();
    }

    // Prepare memory for energy, force
    ViscosityBody visco(3*N,viscosity,f,thermalX,thermalG,prevX);

    // Create the Augmented Lagrangian volume constraint body
    ExactAreaConstraint constraint(N, f, xpos, posGrad, ops.getPolyData());

    // Create Model
    Model model(6*N,f,g);
    model.addBody(&ops);
    model.addBody(&brown);
    model.addBody(&visco);
    model.addBody(&constraint);
    // ****************************************************************//

    // ******************** Read parameter schedule ********************//
    // Input file should contain the following columns
    // Alpha Beta Gamma PercentStrain AreaConstraint NumIterations PrintStep
    sstm << "schedule-" << taskId.str() << ".dat";
    std::ifstream coolFile(sstm.str());
    sstm.str("");
    sstm.clear();
    std::vector<std::vector<double_t> > coolVec;
    double_t currAlpha, currBeta, currGamma, currPercentStrain,
                    currViterMax, currPrintStep, currArea;
    std::string headerline;
    std::getline(coolFile, headerline);// Eat up header line
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

    // ***************** Prepare Output Data files *********************//
    // Identify the Input structure name
    std::string fname = baseFileName;
    std::string dataOutputFile;

    // Detailed output data file
    ofstream detailedOP;
    sstm << fname << "-DetailedOutput-" << taskId.str() << ".dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    detailedOP.open(dataOutputFile.c_str(), std::ofstream::out);
    detailedOP << "#Step" <<"\t"
               << "Beta" << "\t"
               << "Gamma" << "\t"
               << "Asphericity" << "\t"
               << "Radius"  <<"\t"
               << "Volume"  <<"\t"
               << "Area"  <<"\t"
               << "MorseEn"  <<"\t"
               << "NormEn"  <<"\t"
               << "CircEn"  <<"\t"
               << "BrownEn"  <<"\t"
               << "ViscoEn"  <<"\t"
               << "MSD" << "\t"
               << "RMSAngleDeficit"
               << std::endl;

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e5;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, model, f, x, g);
    solver.turnOffLogging();
    // *****************************************************************//

    // ************************ OUTER SOLUTION LOOP **********************//
    // Determine the correct value for the loop start indices
    size_t rowId = 0, colId;
    size_t totSteps = 0;
    while( totSteps < step ){
        totSteps += coolVec[rowId++][5];
    };
    rowId = rowId > 0? rowId - 1 : 0;
    colId = step > 0? coolVec[rowId][5] - (totSteps - step) : 0;

    // The outer loop
    size_t printStep;
    for(auto z=rowId; z < coolVec.size(); ++z){
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
        constraint.setConstraint(constrainedVal);

        // For the very first iteration solve at zero temperature first
        if( z == 0 && colId == 0){
            brown.setCoefficient(0.0);
            visco.setViscosity(0.0);
            solver.solve();
            ops.saveInitialPosition();
            ops.getInitialPositions( initPos );
        }

        // Update prevX
        prevX = x.head(3*N);

        // Set the viscosity and Brownian coefficient
        viscosity = alpha;
        brownCoeff = std::sqrt( 2*alpha/beta );
        brown.setCoefficient(brownCoeff);
        visco.setViscosity(viscosity);

        //**************  INNER SOLUTION LOOP ******************//
        for (auto viter = colId; viter < viterMax; ++viter) {
            // Generate Brownian Kicks
            brown.generateParallelKicks();

            // Set the starting guess for Lambda and K for
            // Augmented Lagrangian
            constraint.setLagrangeCoeff(10.0);
            constraint.setPenaltyCoeff(1000.0);

            // *************** Augmented Lagrangian Loop ************** //
            bool constraintMet = false;
            size_t alIter = 0, alMaxIter = 10;

            while( !constraintMet && (alIter < alMaxIter)){

                // Solve the unconstrained minimization
                solver.solve();

                //Uzawa update
                constraint.uzawaUpdate();

                // Update termination check quantities
                alIter++;
                constraintMet = constraint.constraintSatisfied();
            }
            // *********************************************************//

            // Apply Kabsch Algorithm
            ops.applyKabschAlgorithm();

            // Update kdTree, polyData and neighbors
            ops.updatePolyData();
            ops.updateNeighbors();

            // Check step and save state if needed
            if( step % saveFreq == (saveFreq - 1) ){
                engine = brown.getRandomEngine();
                rng = brown.getRandomGenerator();
                state = SimulationState(N,nameSuffix,step,x,prevX,
                                        initPos,neighbors,engine,rng);
                state.writeToFile(stateFileName);
            }

            // We will print only after every currPrintStep iterations
            if (viter % printStep == 0 && printStep <= viterMax) {
                sstm << fname << "-" << taskId.str() << "-relaxed-"
                     << nameSuffix++ <<".vtk";
                std::string rName = sstm.str();
                ops.printVTKFile(rName);
                sstm.str("");
                sstm.clear();
            }

            std::vector<double_t> msds(2,0);
            msds = ops.getMSD();

            double_t volume = ops.getVolume();
            double_t morseEn = ops.getMorseEnergy();
            double_t normEn = ops.getNormalityEnergy();
            double_t circEn = ops.getCircularityEnergy();

            // Write output to data file
            detailedOP << step << "\t"
                       << beta << "\t"
                       << gamma << "\t"
                       << ops.getAsphericity() << "\t"
                       << ops.getAverageRadius() << "\t"
                       << volume << "\t"
                       << ops.getArea() << "\t"
                       << morseEn << "\t"
                       << normEn << "\t"
                       << circEn << "\t"
                       << brown.getBrownianEnergy() << "\t"
                       << visco.getViscosityEnergy() << "\t"
                       << msds[0] << "\t"
                       << ops.getRMSAngleDeficit()
                       << std::endl;

            // Update prevX
            prevX = x.head(3*N);
            step++;
        }
    }
    // **********************************************************************//

    detailedOP.close();
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 1;
}
