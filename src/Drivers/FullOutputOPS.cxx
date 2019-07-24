#include "ALConstraint.h"
#include "BrownianBody.h"
#include "H5Cpp.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <vtkPolyDataReader.h>

using namespace H5;
using namespace OPS;

#define ARR_SIZE 300

struct DetailedOutput {
  int step;
  int paraviewStep;
  double_t alpha;
  double_t beta;
  double_t gamma;
  double_t asphericity;
  double_t radius;
  double_t volume;
  double_t area;
  double_t morseEn;
  double_t normEn;
  double_t circEn;
  double_t opsEn;
  double_t brownEn;
  double_t viscoEn;
  double_t totEn;
  double_t msd;
  double_t avgDisp;
  double_t maxDisp;
  double_t x[ARR_SIZE];
  double_t y[ARR_SIZE];
  double_t z[ARR_SIZE];
  double_t nx[ARR_SIZE];
  double_t ny[ARR_SIZE];
  double_t nz[ARR_SIZE];
  double_t ndx[ARR_SIZE];
  double_t ndy[ARR_SIZE];
  double_t ndz[ARR_SIZE];

  // Constructor for the struct
  DetailedOutput(int s, int ps, double_t a, double_t b, double_t g, double_t as,
                 double_t r, double_t v, double_t ar, double_t me, double_t ne,
                 double_t ce, double_t oe, double_t be, double_t ve,
                 double_t te, double_t msd, double_t ad, double_t md,
                 Eigen::Matrix3Xd &pos, Eigen::Matrix3Xd &normals,
                 Eigen::Matrix3Xd &normdiff)
      : step(s), paraviewStep(ps), alpha(a), beta(b), gamma(g), asphericity(as),
        radius(r), volume(v), area(ar), morseEn(me), normEn(ne), circEn(ce),
        opsEn(oe), brownEn(be), viscoEn(ve), totEn(te), msd(msd), avgDisp(ad),
        maxDisp(md) {

    for (int i = 0; i < pos.cols(); i++) {
      x[i] = pos(0, i);
      y[i] = pos(1, i);
      z[i] = pos(2, i);
      nx[i] = normals(0, i);
      ny[i] = normals(1, i);
      nz[i] = normals(2, i);
    }
    for (int i = 0; i < normdiff.cols(); i++) {
      ndx[i] = normdiff(0, i);
      ndy[i] = normdiff(1, i);
      ndz[i] = normdiff(2, i);
    }
  }
}; /* ----------  end of struct DetailedOutput  ---------- */

int main(int argc, char *argv[]) {
  clock_t t1, t2, t3;
  t1 = clock();

  if (argc != 2) {
    cout << "usage: " << argv[0] << " <filename>\n";
    return -1;
  }

  bool loggingOn = false;
  bool spikePrintOn = false;

  // ***************** Read Input VTK File *****************//
  std::string inputFileName = argv[1];

  auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
  vtkSmartPointer<vtkPolyData> mesh;

  reader->SetFileName(inputFileName.c_str());
  reader->Update();
  mesh = reader->GetOutput();
  // ********************************************************//

  // ******************* Read Simulation Parameters *********//
  double_t De = 1.0, re = 1.0, s = 7.0;
  double_t alpha = 1.0, beta = 1.0, gamma = 1.0;
  double_t percentStrain = 15;

  std::string constraintType("NULL");
  enum Constraint { AvgArea, AvgVol, ExactArea, ExactVol, ExactAreaAndVolume };
  // int lat_res=100, long_res=101;
  size_t viterMax = 1000;
  size_t nameSuffix = 0;

  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  std::string temp;
  miscInpFile >> temp >> De >> temp >> re >> temp >> constraintType;

  miscInpFile.close();
  s = (100 / (re * percentStrain)) * log(2.0);

  // Validate constraint type
  Constraint type;
  if (constraintType.compare("AverageArea") == 0) {
    type = AvgArea;
  } else if (constraintType.compare("AverageVolume") == 0) {
    type = AvgArea;
  } else if (constraintType.compare("ExactArea") == 0) {
    type = ExactArea;
  } else if (constraintType.compare("ExactVolume") == 0) {
    type = ExactVol;
  } else if (constraintType.compare("ExactAreaAndVolume") == 0) {
    type = ExactAreaAndVolume;
  } else {
    std::cout << "Invalid constraint type specified." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Input file should contain the following columns
  // Alpha Beta Gamma PercentStrain AreaConstraint NumIterations PrintStep
  std::ifstream coolFile("schedule.dat");
  assert(coolFile);
  std::vector<std::vector<double_t>> coolVec;
  double_t currAlpha, currBeta, currGamma, currPercentStrain, currViterMax,
      currPrintStep, currArea;

  std::string headerline;
  std::getline(coolFile, headerline);

  while (coolFile >> currAlpha >> currBeta >> currGamma >> currPercentStrain >>
         currArea >> currViterMax >> currPrintStep) {
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
  // Calculate the number of bonds;
  size_t numBonds = (int)((12 * 5 + (N - 12) * 6) / 2);

  // Generate Rotation Vectors from input point coordinates
  Eigen::Matrix3Xd coords(3, N);
  for (auto i = 0; i < N; ++i) {
    Eigen::Vector3d cp = Eigen::Vector3d::Zero();
    mesh->GetPoint(i, &(cp(0)));
    coords.col(i) = cp;
  }

  // Prepare memory for energy, force
  double_t f;
  Eigen::VectorXd x(6 * N), g(6 * N), prevX(3 * N);
  g.setZero(g.size());
  x.setZero(x.size());

  // Fill x with coords and rotVecs
  Eigen::Map<Eigen::Matrix3Xd> xpos(x.data(), 3, N), xrot(&(x(3 * N)), 3, N),
      prevPos(prevX.data(), 3, N);
  xpos = coords;

  // Renormalize the shell by average edge length
  x /= getPointCloudAvgEdgeLen(inputFileName);
  prevX = x.head(3 * N);
  OPSBody::initialRotationVector(xpos, xrot);
  

  // The average radius
  double_t R0 = xpos.colwise().norm().sum() / N;

  // Create OPSBody
  Eigen::Map<Eigen::Matrix3Xd> posGrad(g.data(), 3, N),
      rotGrad(&g(3 * N), 3, N);
  OPSMesh ops(N, f, R0, xpos, xrot, posGrad, rotGrad);
  ops.setMorseDistance(re);
  ops.setMorseEnergy(De);
  s = 100 * log(2.0) / (re * percentStrain);
  ops.setMorseWellWidth(s);

  // Create Brownian and Viscosity bodies
  Eigen::Map<Eigen::VectorXd> thermalX(x.data(), 3 * N, 1);
  Eigen::Map<Eigen::VectorXd> thermalG(g.data(), 3 * N, 1);
  double_t brownCoeff = 1.0, viscosity = 1.0;
  BrownianBody brown(3 * N, brownCoeff, f, thermalX, thermalG, prevX);
  ViscosityBody visco(3 * N, viscosity, f, thermalX, thermalG, prevX);

  // Create the Augmented Lagrangian volume constraint body
  ALConstraint *constraint;
  if (type == AvgArea) {
    constraint = new AvgAreaConstraint(N, f, xpos, posGrad);
  } else if (type == AvgVol) {
    constraint = new AvgVolConstraint(N, f, xpos, posGrad);
  } else if (type == ExactArea) {
    vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
    constraint = new ExactAreaConstraint(N, f, xpos, posGrad, poly);
  } else if (type == ExactVol) {
    vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
    constraint = new ExactVolConstraint(N, f, xpos, posGrad, poly);
  } else if (type == ExactAreaAndVolume) {
    vtkSmartPointer<vtkPolyData> poly = ops.getPolyData();
    constraint = new ExactAreaVolConstraint(N, f, xpos, posGrad, poly);
    constraint->setTolerance(1e-8);
  }

  // Create Model
  Model model(6 * N, f, g);
  model.addBody(&ops);
  model.addBody(&brown);
  model.addBody(&visco);
  model.addBody(constraint);
  // ****************************************************************//

  // ***************** Prepare Output Data files *********************//
  // Array type for hdf5 compound type
  hsize_t dim[] = {ARR_SIZE};
  ArrayType arr_type(PredType::NATIVE_DOUBLE, 1, dim);

  // Finally create the compound datatype for HDF5
  CompType mtype(sizeof(DetailedOutput));
  mtype.insertMember("Step", HOFFSET(DetailedOutput, step),
                     PredType::NATIVE_INT);
  mtype.insertMember("ParaviewStep", HOFFSET(DetailedOutput, paraviewStep),
                     PredType::NATIVE_INT);
  mtype.insertMember("Alpha", HOFFSET(DetailedOutput, alpha),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Beta", HOFFSET(DetailedOutput, beta),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Gamma", HOFFSET(DetailedOutput, gamma),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Asphericity", HOFFSET(DetailedOutput, asphericity),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Radius", HOFFSET(DetailedOutput, radius),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Volume", HOFFSET(DetailedOutput, volume),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("Area", HOFFSET(DetailedOutput, area),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("MorseEn", HOFFSET(DetailedOutput, morseEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("NormaEn", HOFFSET(DetailedOutput, normEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("CircEn", HOFFSET(DetailedOutput, circEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("OPSEn", HOFFSET(DetailedOutput, opsEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("BrownEn", HOFFSET(DetailedOutput, brownEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("ViscoEn", HOFFSET(DetailedOutput, viscoEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("TotalEn", HOFFSET(DetailedOutput, totEn),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("MSD", HOFFSET(DetailedOutput, msd),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("AvgDisp", HOFFSET(DetailedOutput, avgDisp),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("MaxDisp", HOFFSET(DetailedOutput, maxDisp),
                     PredType::NATIVE_DOUBLE);
  mtype.insertMember("X", HOFFSET(DetailedOutput, x), arr_type);
  mtype.insertMember("Y", HOFFSET(DetailedOutput, y), arr_type);
  mtype.insertMember("Z", HOFFSET(DetailedOutput, z), arr_type);
  mtype.insertMember("NX", HOFFSET(DetailedOutput, nx), arr_type);
  mtype.insertMember("NY", HOFFSET(DetailedOutput, ny), arr_type);
  mtype.insertMember("NZ", HOFFSET(DetailedOutput, nz), arr_type);
  mtype.insertMember("NDX", HOFFSET(DetailedOutput, ndx), arr_type);
  mtype.insertMember("NDY", HOFFSET(DetailedOutput, ndy), arr_type);
  mtype.insertMember("NDZ", HOFFSET(DetailedOutput, ndz), arr_type);

  // Create a dataspace for our dataset with Compound datatype
  hsize_t startDim[] = {1};
  hsize_t maxDim[] = {H5S_UNLIMITED};
  DataSpace memspace(1, startDim, maxDim);

  // Create property list to enable chunking and compression in our dataset
  hsize_t chunk[] = {10000};
  DSetCreatPropList dsprop;
  dsprop.setChunk(1, chunk);
  // dsprop.setDeflate(4);

  // Identify the Input structure name
  std::string fname = inputFileName.substr(0, inputFileName.find("."));
  std::stringstream sstm;
  std::string dataOutputFile;
  sstm << fname << "-DetailedOutput.h5";
  dataOutputFile = sstm.str();
  sstm.str("");
  sstm.clear();

  // Create a file to create a dataset into
  H5File file(dataOutputFile, H5F_ACC_TRUNC);

  // Create a dataset in the file using above dataspace and datatype
  DataSet data = file.createDataSet("DetailedOutput", mtype, memspace, dsprop);

  // Create attributes that contain the numBonds and numPoints
  H5std_string P("NumParticles");
  H5std_string B("NumBonds");
  hsize_t attrDim[] = {1};
  DataSpace attrDS(1, attrDim);
  Attribute particles = data.createAttribute(P, PredType::NATIVE_UINT, attrDS);
  Attribute bonds = data.createAttribute(B, PredType::NATIVE_UINT, attrDS);
  particles.write(PredType::NATIVE_UINT, &N);
  bonds.write(PredType::NATIVE_UINT, &numBonds);

  // Create the output file for average data
  ofstream outerLoopFile;
  sstm << fname << "-AverageOutput.dat";
  dataOutputFile = sstm.str();
  sstm.str("");
  sstm.clear();
  outerLoopFile.open(dataOutputFile.c_str());
  outerLoopFile << "#BigStep"
                << "\t"
                << "PercentStrain"
                << "\t"
                << "Alpha"
                << "\t"
                << "Beta"
                << "\t"
                << "Gamma"
                << "\t"
                << "PercentStrain"
                << "\t"
                << "Radius"
                << "\t"
                << "Asphericity" << std::endl;
  // ******************************************************************//

  // ************************* Create Solver ************************  //
  size_t m = 5, iprint = 1000, maxIter = 1e5;
  double_t factr = 10.0, pgtol = 1e-8;
  LBFGSBParams solverParams(m, iprint, maxIter, factr, pgtol);
  LBFGSBWrapper solver(solverParams, model, f, x, g);
  solver.turnOffLogging();
  // *****************************************************************//

  // ********************* Prepare data for simulation ****************//
  // Calculate Average Edge Length
  double_t avgEdgeLen = ops.getAverageEdgeLength();
  if (loggingOn) {
    std::cout << "Initial Avg Edge Length = " << avgEdgeLen << std::endl;
  }
  // Renormalize positions such that avgEdgeLen = 1.0
  for (auto i = 0; i < N; ++i) {
    xpos.col(i) = xpos.col(i) / avgEdgeLen;
  }

  // Update the OPSBody member variables as per new positions
  ops.updatePolyData();
  ops.updateNeighbors();
  ops.saveInitialPosition(); /*!< For Mean Squared Displacement */
  avgEdgeLen = ops.getAverageEdgeLength();
  if (loggingOn)
    std::cout << "After renormalizing, Avg Edge Length = " << avgEdgeLen
              << std::endl;
  // ******************************************************************//

  t3 = clock();
  // ************************ OUTER SOLUTION LOOP **********************//
  int printStep;
  hsize_t step = 0;
  int paraviewStep = -1;
  for (int z = 0; z < coolVec.size(); z++) {
    alpha = coolVec[z][0];
    beta = coolVec[z][1];
    gamma = coolVec[z][2];
    percentStrain = coolVec[z][3];
    double_t constrainedVal = coolVec[z][4];
    viterMax = coolVec[z][5];
    printStep = (int)coolVec[z][6];

    // Update OPS params
    s = (100 / (avgEdgeLen * percentStrain)) * log(2.0);
    ops.setFVK(gamma);
    ops.setMorseWellWidth(s);

    // Set up the constraint value as the zero temperature value
    constraint->setConstraint(constrainedVal);

    // For the very first iteration solve at zero temperature first
    if (z == 0) {
      brown.setCoefficient(0.0);
      visco.setViscosity(0.0);
      solver.solve();
    }

    // Update prevX
    prevX = x.head(3 * N);

    // Set the viscosity and Brownian coefficient
    viscosity = alpha * De / (avgEdgeLen * avgEdgeLen);
    brownCoeff = std::sqrt(2 * alpha / beta) * (De / avgEdgeLen);
    if (loggingOn) {
      std::cout << "Viscosity = " << viscosity << std::endl;
      std::cout << "Brownian Coefficient = " << brownCoeff << std::endl;
    }
    brown.setCoefficient(brownCoeff);
    visco.setViscosity(viscosity);

    //**************  INNER SOLUTION LOOP ******************//
    Eigen::Matrix3Xd averagePosition(3, N);
    averagePosition = Eigen::Matrix3Xd::Zero(3, N);

    // Average energy across time steps
    double_t avgTotalEnergy = 0.0;

    for (int viter = 0; viter < viterMax; viter++) {
      if (loggingOn)
        std::cout << std::endl
                  << "VISCOUS ITERATION: " << step << std::endl
                  << std::endl;

      // Generate Brownian Kicks
      brown.generateParallelKicks();

      // Store data for Kabsch
      ops.updateDataForKabsch();

      // Set the starting guess for Lambda and K for
      // Augmented Lagrangian
      constraint->setLagrangeCoeff(10.0);
      constraint->setPenaltyCoeff(1000.0);

      // *************** Augmented Lagrangian Loop ************** //
      bool constraintMet = false;
      size_t alIter = 0, alMaxIter = 10;

      while (!constraintMet && (alIter < alMaxIter)) {
        if (loggingOn)
          std::cout << "Augmented Lagrangian iteration: " << alIter
                    << std::endl;

        // Solve the unconstrained minimization
        solver.solve();

        // Uzawa update
        constraint->uzawaUpdate();

        // Update termination check quantities
        alIter++;
        constraintMet = constraint->constraintSatisfied();
      }
      if (loggingOn) {
        constraint->printCompletion();
        std::cout << "Constraint satisfied in " << alIter << " iterations."
                  << std::endl
                  << std::endl;
      }
      // *********************************************************//

      // Apply Kabsch Algorithm
      ops.applyKabschAlgorithm();

      // Update kdTree, polyData and neighbors
      ops.updatePolyData();
      ops.updateNeighbors();

      // Add current solution to average position data
      averagePosition += xpos;

      // Calculate statistics about average displacement and
      // largest displacement
      double_t avgDisplacement, maxDisplacement;
      Eigen::Matrix3Xd posDiff = xpos - prevPos;
      avgDisplacement = posDiff.colwise().norm().sum() / N;
      maxDisplacement = posDiff.colwise().norm().maxCoeff();
      Eigen::Matrix3Xd normDiff(3, numBonds);
      ops.getDiffNormals(normDiff);
      Eigen::Matrix3Xd normals(3, numBonds);
      ops.getNormals(normals);

      //********** Print relaxed configuration ************//
      // We will print only after every currPrintStep iterations
      if (viter % printStep == 0 && printStep < viterMax) {
        paraviewStep++;
        sstm << fname << "-relaxed-" << nameSuffix++ << ".vtk";
        std::string rName = sstm.str();
        ops.printVTKFile(rName);
        sstm.str("");
        sstm.clear();
      }
      // Print VTK file if there is an abrupt change in energy
      if (spikePrintOn) {
        if (std::abs(avgTotalEnergy) > 0 &&
            std::abs((f - avgTotalEnergy) / avgTotalEnergy) > 3) {
          sstm << fname << "-Spike-" << step << ".vtk";
          std::string rName = sstm.str();
          ops.printVTKFile(rName);
          sstm.str("");
          sstm.clear();
        }
      }

      int paraviewStepPrint;
      paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;

      DetailedOutput dop(step, paraviewStepPrint, alpha, beta, gamma,
                         ops.getAsphericity(), ops.getAverageRadius(),
                         ops.getVolume(), ops.getArea(), ops.getMorseEnergy(),
                         ops.getNormalityEnergy(), ops.getCircularityEnergy(),
                         ops.getTotalEnergy(), brown.getBrownianEnergy(),
                         visco.getViscosityEnergy(), f,
                         ops.getMeanSquaredDisplacement(), avgDisplacement,
                         maxDisplacement, posDiff, normals, normDiff);

      // Write data to the HDF5 file
      hsize_t size[] = {step + 1};
      hsize_t start[] = {step};
      hsize_t count[] = {1};
      data.extend(size);
      DataSpace fspace = data.getSpace();
      fspace.selectHyperslab(H5S_SELECT_SET, count, start);
      data.write(&dop, mtype, memspace, fspace);

      // Update prevX
      prevX = x.head(3 * N);
      avgTotalEnergy = (avgTotalEnergy * step + f) / (step + 1);
      std::cout << " Average Total Energy = " << avgTotalEnergy << std::endl;
      step++;
    }
    //************************************************//

    // Calculate the average particle positions and avg radius
    double avgShapeRad = 0.0, avgShapeAsph = 0.0;
    auto avgPos = vtkSmartPointer<vtkPoints>::New();
    auto avgPosData = vtkSmartPointer<vtkDoubleArray>::New();
    averagePosition = averagePosition / viterMax;
    avgShapeRad = averagePosition.colwise().norm().sum() / N;
    void *avgPosPtr = (void *)averagePosition.data();
    avgPosData->SetVoidArray(avgPosPtr, 3 * N, 1);
    avgPosData->SetNumberOfComponents(3);
    avgPos->SetData(avgPosData);

    sstm << fname << "-AvgShape-" << z << ".vtk";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();

    // Print the average shape
    delaunay3DSurf(avgPos, dataOutputFile);

    // Calculate the asphericity of the average shape
    Eigen::RowVectorXd R(N);
    R = averagePosition.colwise().norm();
    avgShapeAsph += ((R.array() - avgShapeRad).square()).sum();
    avgShapeAsph /= (N * avgShapeRad * avgShapeRad);

    outerLoopFile << z << "\t" << percentStrain << "\t" << alpha << "\t" << beta
                  << "\t" << gamma << "\t" << percentStrain << "\t"
                  << avgShapeRad << "\t" << avgShapeAsph << std::endl;
  }
  // *****************************************************************************//

  file.close();
  outerLoopFile.close();
  t2 = clock();
  float diff((float)t2 - (float)t1);
  std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
            << " seconds" << std::endl;
  delete constraint;
  return 1;
}
