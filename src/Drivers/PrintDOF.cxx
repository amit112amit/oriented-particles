#include "HelperFunctions.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSBody.h"
#include "OPSModel.h"
#include <random>
#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>

using namespace OPS;
using namespace std;

int main(int argc, char *argv[]) {
  clock_t t1;
  t1 = clock();

  if (argc != 2) {
    cout << "usage: " << argv[0] << " <filename>\n";
    return -1;
  }

  //******************** Create variables ********************//
  std::string inFile = argv[1];
  std::string baseFileName = inFile.substr(0, inFile.length() - 4);

  double_t alpha = 1.0, beta = 1.0, gamma = 1.0, re = 1.0, f = 0, R0 = 0,
           percentStrain = 15, s = (100 / (re * percentStrain)) * log(2.0);

  size_t viterMax, nameSuffix = 0, step = 0, N = 10, saveFreq = 100000;

  Eigen::VectorXd x(6 * N), prevX(3 * N), g(6 * N);
  Eigen::Matrix3Xd initPos(3, N);
  std::vector<size_t> neighbors(N);
  std::mt19937 engine;
  auto rng = std::normal_distribution<double_t>(0.0, 1.0);
  //**********************************************************//

  //******************** Resume or fresh start? ********************//
  // Create a SimulationState
  SimulationState state;

  // Check if a state file exists in current directory
  std::string stateFileName("SimulationState.dat");
  ifstream stateFile(stateFileName.c_str());

  if (stateFile.good()) {
    // Read previously stored state to resume the simulation
    state = SimulationState::readFromFile(stateFileName);
    nameSuffix = state.getNameSuffix();
    nameSuffix = nameSuffix > 0 ? nameSuffix + 1 : 0;
    step = state.getStep() + 1;
    N = state.getN();
    engine = state.getRandomEngine();
    rng = state.getRandomGenerator();
    R0 = state.getRadius0();
    // Resize matrices and vectors
    x.resize(6 * N);
    g.resize(6 * N);
    prevX.resize(3 * N);
    initPos.resize(3, N);
    neighbors.resize(N);
    // Read the vectors into local variables
    x = state.getX();
    prevX = state.getPrevX();
    neighbors = state.getNeighbors();
    initPos = state.getInitPos();
  } else {
    // No state file found. So we need to start a new simulation
    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(inFile.c_str());
    reader->ReadAllVectorsOn();
    reader->Update();
    auto mesh = reader->GetOutput();
    N = mesh->GetNumberOfPoints();
    // Resize matrices and vectors
    x.resize(6 * N);
    g.resize(6 * N);
    prevX.resize(3 * N);
    neighbors.resize(N);
    // Read point coordinates from input mesh
    for (auto i = 0; i < N; ++i)
      mesh->GetPoint(i, &x(3 * i));
    // Renormalize by the average edge length
    x /= getPointCloudAvgEdgeLen(inFile);

    Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(), 3, N), xrot(&x(3 * N), 3, N);

    // The average radius
    R0 = xpos.colwise().norm().sum() / N;

    // Generate rotation vectors from input point coordinates
    OPSBody::initialRotationVector(xpos, xrot);
  }
  // ****************************************************************//

  // ***************** Create Bodies and Model ****************//

  // Create OPSModel
  g.setZero(g.size());
  OPSModel ops(N, f, x, g, prevX, R0);
  ops.setMorseDistance(re);
  ops.setMorseWellWidth(s);
  ops.updateTriangles();

  // Set ops history variables
  if (stateFile.good()) {
    ops.setInitialNeighbors(neighbors);
    ops.setInitialPositions(initPos);
    ops.setRandomEngine(engine);
    ops.setRandomGenerator(rng);
  } else {
    prevX = x.head(3 * N);
    neighbors = ops.getInitialNeighbors();
  }

  // Create Model
  Model model(6 * N, f, g);
  model.addBody(&ops);
  // ****************************************************************//

  // ******************** Read parameter schedule ********************//
  // Input file should contain the following columns
  // Alpha Beta Gamma PercentStrain AreaConstraint NumIterations PrintStep
  std::ifstream coolFile("schedule.dat");
  assert(coolFile);
  std::vector<std::vector<double_t>> coolVec;
  double_t currAlpha, currBeta, currGamma, currPercentStrain, currViterMax,
      currArea;
  std::string headerline;
  std::getline(coolFile, headerline); // Eat up header line
  while (coolFile >> currAlpha >> currBeta >> currGamma >> currPercentStrain >>
         currArea >> currViterMax) {
    std::vector<double> currLine;
    currLine.push_back(currAlpha);
    currLine.push_back(currBeta);
    currLine.push_back(currGamma);
    currLine.push_back(currPercentStrain);
    currLine.push_back(currArea);
    currLine.push_back(currViterMax);
    coolVec.push_back(currLine);
  }
  coolFile.close();
  // **********************************************************//

  // ***************** Prepare Output Data files *********************//
  // Identify the Input structure name
  std::string fname = baseFileName;
  std::stringstream sstm;
  std::string dataOutputFile;

  // Detailed output data file
  ofstream detailedOP;
  sstm << fname << "-DOF.dat";
  dataOutputFile = sstm.str();
  sstm.str("");
  sstm.clear();
  detailedOP.open(dataOutputFile);
  detailedOP.precision(5);

  // ************************* Create Solver ************************  //
  size_t m = 5, iprint = 1000, maxIter = 1e5;
  double_t factr = 10.0, pgtol = 1e-8;
  LBFGSBParams solverParams(m, iprint, maxIter, factr, pgtol);
  LBFGSBWrapper solver(solverParams, model, f, x, g);
  solver.turnOffLogging();
  // *****************************************************************//

  // ************************ OUTER SOLUTION LOOP **********************//
  // Determine the correct value for the loop start indices
  size_t rowId = 0, colId;
  size_t totSteps = 0;
  while (totSteps < step) {
    totSteps += coolVec[rowId++][5];
  };
  rowId = rowId > 0 ? rowId - 1 : 0;
  colId = step > 0 ? coolVec[rowId][5] - (totSteps - step) : 0;

  // The outer loop
  for (auto z = rowId; z < coolVec.size(); ++z) {
    alpha = coolVec[z][0];
    beta = coolVec[z][1];
    gamma = coolVec[z][2];
    percentStrain = coolVec[z][3];
    double_t constrainedVal = coolVec[z][4];
    viterMax = coolVec[z][5];

    // Write gamma and beta to output file
    detailedOP << "#N " << N << " Gamma " << gamma << " Beta " << beta
               << std::endl;

    // Update OPS params
    s = (100 / percentStrain) * log(2.0);
    ops.setFVK(gamma);
    ops.setMorseWellWidth(s);

    // Set up the constraint value as the zero temperature value
    ops.setConstraint(constrainedVal);

    // For the very first iteration solve at zero temperature first
    if (z == 0 && colId == 0) {
      ops.setBrownCoeff(0.0);
      ops.setViscosity(0.0);
      solver.solve();
      ops.saveInitialPosition();
      ops.getInitialPositions(initPos);
      // Write the zero temperature x to data file
      for (auto idx = 0; idx < 6 * N - 1; ++idx)
        detailedOP << x(idx) << ",";
      detailedOP << x(6 * N - 1) << std::endl;
    }

    // Update prevX
    prevX = x.head(3 * N);

    // Set the viscosity and Brownian coefficient
    double_t viscosity = alpha;
    double_t brownCoeff = std::sqrt(2 * alpha / beta);
    ops.setBrownCoeff(brownCoeff);
    ops.setViscosity(viscosity);

    //**************  INNER SOLUTION LOOP ******************//
    for (auto viter = colId; viter < viterMax; ++viter) {
      // Generate Brownian Kicks
      ops.generateParallelKicks();

      // Set the starting guess for Lambda and K for
      // Augmented Lagrangian
      ops.setLagrangeCoeff(10.0);
      ops.setPenaltyCoeff(1000.0);

      // *************** Augmented Lagrangian Loop ************** //
      bool constraintMet = false;
      size_t alIter = 0, alMaxIter = 10;

      while (!constraintMet && (alIter < alMaxIter)) {

        // Solve the unconstrained minimization
        solver.solve();

        // Uzawa update
        ops.uzawaUpdate();

        // Update termination check quantities
        alIter++;
        constraintMet = ops.constraintSatisfied();
      }
      // *********************************************************//

      // Apply Kabsch Algorithm
      ops.applyKabschAlgorithm();

      // Update kdTree, polyData and neighbors
      ops.updateTriangles();

      // Check step and save state if needed
      if (step % saveFreq == (saveFreq - 1)) {
        engine = ops.getRandomEngine();
        rng = ops.getRandomGenerator();
        state = SimulationState(N, nameSuffix, step, gamma, R0, beta, x, prevX,
                                initPos, neighbors, engine, rng);
        state.writeToFile("SimulationState.dat");
      }

      // Write x to data file
      for (auto idx = 0; idx < 6 * N - 1; ++idx)
        detailedOP << x(idx) << ",";
      detailedOP << x(6 * N - 1) << std::endl;

      // Update prevX
      prevX = x.head(3 * N);
      step++;
    }
  }
  // **********************************************************************//

  detailedOP.close();
  float diff((float)clock() - (float)t1);
  std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC << " seconds"
            << std::endl;
  return 1;
}
