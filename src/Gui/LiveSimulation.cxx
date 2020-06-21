#include "LiveSimulation.h"

namespace OPS
{

void LiveSimulation::Initialize()
{
  // Initialize circular buffer with capacity 50000
  _rmsAD = new QwtCircBuffSeriesData(5000);
  _opsEn = new QwtCircBuffSeriesData(5000);
  _vol = new QwtCircBuffSeriesData(5000);

  // ***************** Read Input VTK File *****************//
  std::string inputFileName = _inputVTK;

  auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
  vtkSmartPointer<vtkPolyData> mesh;

  reader->SetFileName(inputFileName.c_str());
  reader->ReadAllVectorsOn();
  reader->Update();
  mesh = reader->GetOutput();
  // ********************************************************//

  // ******************* Simulation Parameters *********//
  double_t re = 1.0;
  double_t percentStrain = 15;
  double_t s = (100 / (re * percentStrain)) * log(2.0);
  double_t brownCoeff = 1.0, viscosity = 1.0;
  ReadFvKAreaData(_inputZeroData);
  _zeroOpsEnVal = GetInterpolatedValue(_gamma, _opsEnDat);
  _zeroRmsAdVal = GetInterpolatedValue(_gamma, _rmsAdDat);
  _zeroVolVal = GetInterpolatedValue(_gamma, _volDat);

  // **********************************************************//

  // ***************** Create Bodies and Model ****************//
  // Set number of OPS particles
  _N = mesh->GetNumberOfPoints();

  // Read point coordinates from input mesh
  Eigen::Matrix3Xd coords(3, _N);
  for (auto i = 0; i < _N; ++i)
  {
    Eigen::Vector3d cp = Eigen::Vector3d::Zero();
    mesh->GetPoint(i, &(cp(0)));
    coords.col(i) = cp;
  }

  // Prepare memory for energy, force
  _x = Eigen::VectorXd(6 * _N);
  _g = Eigen::VectorXd(6 * _N);
  _prevX = Eigen::VectorXd(3 * _N);
  _initialX = Eigen::VectorXd(6 * _N);
  _x.setZero(_x.size());
  _g.setZero(_g.size());
  _prevX.setZero(_prevX.size());
  _initialX.setZero(_initialX.size());

  // Fill _x with coords and rotVecs
  Eigen::Map<Eigen::Matrix3Xd> xpos(_x.data(), 3, _N), prevPos(_prevX.data(), 3, _N);
  xpos = coords;

  // Renormalize by the average edge length
  _x /= getPointCloudAvgEdgeLen(inputFileName);
  _prevX = _x.head(3 * _N);

  // Generate rotation vectors from input point coordinates
  Eigen::Map<Eigen::Matrix3Xd> xrot(&(_x(3 * _N)), 3, _N);
  OPSBody::initialRotationVector(xpos, xrot);

  // Calculate the starting average radius
  _R0 = xpos.colwise().norm().sum() / _N;

  // Create OPSBody
  Eigen::Map<Eigen::Matrix3Xd> posGrad(_g.data(), 3, _N),
      rotGrad(&_g(3 * _N), 3, _N);
  _ops = new OPSMesh(_N, _f, _R0, xpos, xrot, posGrad, rotGrad, prevPos);
  _ops->setFVK(_gamma);
  _ops->setMorseDistance(re);
  _ops->setMorseWellWidth(s);

  // Create InternalPressure body
  _pressureBody =
      new PressureBody(_N, _f, xpos, prevPos, posGrad, _ops->getPolyData());

  // Create Brownian and Viscosity bodies
  Eigen::Map<Eigen::VectorXd> thermalX(_x.data(), 3 * _N, 1);
  Eigen::Map<Eigen::VectorXd> thermalG(_g.data(), 3 * _N, 1);
  _brown = new BrownianBody(3 * _N, brownCoeff, _f, thermalX, thermalG, _prevX);
  _visco = new ViscosityBody(3 * _N, viscosity, _f, thermalX, thermalG, _prevX);

  // Create area constraint
  vtkSmartPointer<vtkPolyData> poly = _ops->getPolyData();
  _constraint = new ExactAreaConstraint(_N, _f, xpos, posGrad, poly);
  _constraint->setConstraint(GetInterpolatedValue(_gamma, _areaDat));

  // Create Model
  _model = new Model(6 * _N, _f, _g);
  _model->addBody(_ops);
  _model->addBody(_brown);
  _model->addBody(_visco);
  _model->addBody(_constraint);
  _model->addBody(_pressureBody);
  // ****************************************************************//

  // ***************** Prepare Output Data files *********************//
  // Identify the Input structure name
  std::string fname = inputFileName.substr(0, inputFileName.length() - 4);
  std::stringstream sstm;

  // Detailed output data file
  sstm << fname << "-DetailedOutput.dat";
  _dataOutputFile = sstm.str();
  sstm.str("");
  sstm.clear();
  _detailedOP.open(_dataOutputFile.c_str(), std::ofstream::out);
  _detailedOP << "Alpha"
              << "\t"
              << "Beta"
              << "\t"
              << "Gamma"
              << "\t"
              << "Asphericity"
              << "\t"
              << "MorseEn"
              << "\t"
              << "NormalityEn"
              << "\t"
              << "CircularityEn"
              << "\t"
              << "BrownEn"
              << "\t"
              << "ViscoEn"
              << "\t"
              << "TotalEnergy"
              << "\t"
              << "PressureWork"
              << "\t"
              << "MSD"
              << "\t"
              << "RMSAngleDeficit" << std::endl;

  // ************************* Create Solver ************************  //
  size_t m = 5, iprint = 1000, maxIter = 1e5;
  double_t factr = 10.0, pgtol = 1e-8;
  LBFGSBParams solverParams(m, iprint, maxIter, factr, pgtol);
  _solver = new LBFGSBWrapper(solverParams, *(_model), _f, _x, _g);
  _solver->turnOffLogging();
  // *****************************************************************//

  // ********************* Prepare data for simulation ****************//
  // Update the OPSBody member variables as per new positions
  _ops->updatePolyData();
  _ops->updateNeighbors();
  _ops->saveInitialPosition();
  double_t avgEdgeLen = _ops->getAverageEdgeLength();

  // Relax at zero temperature once
  _brown->setCoefficient(0.0);
  _visco->setViscosity(0.0);
  _constraint->setLagrangeCoeff(0.0);
  _constraint->setPenaltyCoeff(0.0);
  _solver->solve();
  _prevX = _x.head(3 * _N);

  // Save the current _x for resetting the simulation
  _initialX = _x;

  // Set Finite temperature coefficients
  viscosity = _alpha / (avgEdgeLen * avgEdgeLen);
  brownCoeff = std::sqrt(2 * _alpha / _beta) / avgEdgeLen;
  _brown->setCoefficient(brownCoeff);
  _visco->setViscosity(viscosity);
  // ******************************************************************//
  emit simulationReady();
}

//! Load simulation state from a SimulationState object
void LiveSimulation::LoadState(QString st)
{
  // Delete existing objects assigned using new
  delete _ops, _brown, _visco, _constraint, _pressureBody, _model, _solver;

  // Read the state file
  SimulationState state = SimulationState::readFromFile(st.toStdString());
  _N = state.getN();
  _step = state.getStep() + 1;
  _resetStepVal = _step;
  _gamma = state.getGamma();
  _beta = state.getBeta();
  _R0 = state.getRadius0();
  _alpha = 2.5e5;
  std::vector<size_t> neighbors(_N);
  neighbors = state.getNeighbors();

  // Reset the circular buffers
  _rmsAD->clear();
  _opsEn->clear();
  _vol->clear();

  // ******************* Simulation Parameters *********//
  double_t re = 1.0;
  double_t percentStrain = 15;
  double_t s = (100 / (re * percentStrain)) * log(2.0);
  double_t brownCoeff = std::sqrt(2 * _alpha / _beta);
  double_t viscosity = _alpha;
  _zeroOpsEnVal = GetInterpolatedValue(_gamma, _opsEnDat);
  _zeroRmsAdVal = GetInterpolatedValue(_gamma, _rmsAdDat);
  _zeroVolVal = GetInterpolatedValue(_gamma, _volDat);
  // **********************************************************//

  // ***************** Create Bodies and Model ****************//
  // Resize matrices and vectors
  _x = Eigen::VectorXd(6 * _N);
  _g = Eigen::VectorXd(6 * _N);
  _prevX = Eigen::VectorXd(3 * _N);
  _initialX = Eigen::VectorXd(6 * _N);
  _x = state.getX();
  _prevX = state.getPrevX();
  _initialX = _x;

  // Fill coordinates and rotation vectors and prepare gradient vectors
  Eigen::Map<Eigen::Matrix3Xd> xpos(_x.data(), 3, _N),
      xrot(&(_x(3 * _N)), 3, _N), prevPos(_prevX.data(), 3, _N);
  Eigen::Map<Eigen::Matrix3Xd> posGrad(_g.data(), 3, _N),
      rotGrad(&_g(3 * _N), 3, _N);

  // Create OPSBody
  _ops = new OPSMesh(_N, _f, _R0, xpos, xrot, posGrad, rotGrad, prevPos);
  _ops->setFVK(_gamma);
  _ops->setMorseDistance(re);
  _ops->setMorseWellWidth(s);
  _ops->setInitialNeighbors(neighbors);
  _ops->setInitialPositions(state.getInitPos());
  _ops->updateNeighbors();
  vtkSmartPointer<vtkPolyData> poly = _ops->getPolyData();

  // Create InternalPressure body
  _pressureBody = new PressureBody(_N, _f, xpos, prevPos, posGrad, poly);

  // Create Brownian body
  Eigen::Map<Eigen::VectorXd> thermalX(_x.data(), 3 * _N, 1);
  Eigen::Map<Eigen::VectorXd> thermalG(_g.data(), 3 * _N, 1);
  _brown = new BrownianBody(3 * _N, brownCoeff, _f, thermalX, thermalG, _prevX);
  _brown->setRandomEngine(state.getRandomEngine());
  _brown->setRandomGenerator(state.getRandomGenerator());

  // Create viscosity body
  _visco = new ViscosityBody(3 * _N, viscosity, _f, thermalX, thermalG, _prevX);

  // Create area constraint
  _constraint = new ExactAreaConstraint(_N, _f, xpos, posGrad, poly);
  _constraint->setConstraint(GetInterpolatedValue(_gamma, _areaDat));

  // Create Model
  _model = new Model(6 * _N, _f, _g);
  _model->addBody(_ops);
  _model->addBody(_brown);
  _model->addBody(_visco);
  _model->addBody(_constraint);
  _model->addBody(_pressureBody);
  // ****************************************************************//

  // ***************** Prepare Output Data files *********************//
  // Identify the Input structure name
  std::string fname = _inputVTK.substr(0, _inputVTK.length() - 4);
  std::stringstream sstm;

  // Detailed output data file
  sstm << fname << "-DetailedOutput.dat";
  _dataOutputFile = sstm.str();
  sstm.str("");
  sstm.clear();
  _detailedOP.open(_dataOutputFile.c_str(), std::ofstream::out);
  _detailedOP << "Alpha"
              << "\t"
              << "Beta"
              << "\t"
              << "Gamma"
              << "\t"
              << "Asphericity"
              << "\t"
              //<< "MorseEn"
              //<< "\t"
              //<< "NormalityEn"
              //<< "\t"
              //<< "CircularityEn"
              //<< "\t"
              //<< "BrownEn"
              //<< "\t"
              //<< "ViscoEn"
              //<< "\t"
              << "TotalEnergy"
              << "\t"
              << "PressureWork"
              << "\t"
              << "MSD"
              << "\t"
              << "RMSAngleDeficit" << std::endl;

  // ************************* Create Solver ************************  //
  size_t m = 5, iprint = 1000, maxIter = 1e5;
  double_t factr = 10.0, pgtol = 1e-8;
  LBFGSBParams solverParams(m, iprint, maxIter, factr, pgtol);
  _solver = new LBFGSBWrapper(solverParams, *(_model), _f, _x, _g);
  _solver->turnOffLogging();
  // *****************************************************************//

  // Send signals to update the GUI
  emit updatePlotXAxis(_step);
  emit simulationReady();
}

//! Save Simulation State
void LiveSimulation::SaveState(QString s)
{
  Engine engine = _brown->getRandomEngine();
  NormD rng = _brown->getRandomGenerator();
  Eigen::Matrix3Xd initPos(3, _N);
  _ops->getInitialPositions(initPos);
  std::vector<size_t> neighbors = _ops->getInitialNeighbors();
  SimulationState state = SimulationState(
      _N, 0, _step, _gamma, _beta, _R0, _x, _prevX, initPos, neighbors, engine, rng);
  state.writeToFile(s.toStdString());
}

//! Interpolation of area
double_t LiveSimulation::GetInterpolatedValue(double_t x,
                                              std::vector<double_t> &y)
{
  int size = _gammaDat.size();
  int i = 0;
  if (x >= _gammaDat[size - 2])
  {
    i = size - 2;
  }
  else
  {
    while (x > _gammaDat[i + 1])
      i++;
  }
  double_t xL = _gammaDat[i], yL = y[i], xR = _gammaDat[i + 1], yR = y[i + 1];
  if (x < xL)
    yR = yL;
  if (x > xR)
    yL = yR;
  double_t dydx = (yR - yL) / (xR - xL);
  return yL + dydx * (x - xL);
}

// Read gamma and area data from file
void LiveSimulation::ReadFvKAreaData(std::string s)
{
  std::ifstream inputFile(s.c_str());
  std::string line;
  double_t gamma, area, volume, rmsAd, opsEn, ignore;
  // Clean the existing vectors
  _gammaDat.clear();
  _volDat.clear();
  _areaDat.clear();
  _opsEnDat.clear();
  _rmsAdDat.clear();
  // Eat up the header line
  std::getline(inputFile, line);
  while (std::getline(inputFile, line))
  {
    std::istringstream ss(line);
    ss >> gamma >> ignore >> ignore >> volume >> area >> ignore >> ignore >>
        ignore >> opsEn >> rmsAd;
    _gammaDat.push_back(gamma);
    _volDat.push_back(volume);
    _areaDat.push_back(area);
    _opsEnDat.push_back(opsEn);
    _rmsAdDat.push_back(rmsAd);
  }
  // Reverse the vectors so that gamma values are increasing
  std::reverse(_gammaDat.begin(), _gammaDat.end());
  std::reverse(_volDat.begin(), _volDat.end());
  std::reverse(_areaDat.begin(), _areaDat.end());
  std::reverse(_opsEnDat.begin(), _opsEnDat.end());
  std::reverse(_rmsAdDat.begin(), _rmsAdDat.end());
}

// Update gamma and area constraint value
void LiveSimulation::UpdateGamma(double g)
{
  _gamma = g;
  _ops->setFVK(_gamma);
  _constraint->setConstraint(GetInterpolatedValue(_gamma, _areaDat));
  emit updateZeroOpsEn(GetInterpolatedValue(_gamma, _opsEnDat));
  emit updateZeroRmsAd(GetInterpolatedValue(_gamma, _rmsAdDat));
  emit updateZeroVolume(GetInterpolatedValue(_gamma, _volDat));
}

// Update pressure value
void LiveSimulation::UpdatePressure(double p)
{
  _pressure = p;
  _pressureBody->setPressure(p);
}

// Update beta and viscosity and brownian coefficients
void LiveSimulation::UpdateBeta(double b)
{
  if (b < 1e-6)
  {
    _alpha = 0;
    _constraint->setLagrangeCoeff(0.0);
    _constraint->setPenaltyCoeff(0.0);
  }
  else
  {
    _alpha = 2.5e5;
    _beta = 1.0 / b;
  }
  _brown->setCoefficient(std::sqrt(2 * _alpha * b));
  _visco->setViscosity(_alpha);
}

// Start running
void LiveSimulation::SolveOneStep()
{
  if (_keepRunning)
  {
    if (_step == 0)
    {
      // Relax at zero temperature once
      _brown->setCoefficient(0.0);
      _visco->setViscosity(0.0);
      _pressureBody->setPressure(0.0);
      _constraint->setLagrangeCoeff(0.0);
      _constraint->setPenaltyCoeff(0.0);
      _solver->solve();
      _prevX = _x.head(3 * _N);
      _ops->updatePolyData();
      _ops->updateNeighbors();
      _brown->setCoefficient(std::sqrt(2 * _alpha / _beta));
      _visco->setViscosity(_alpha);
      _pressureBody->setPressure(_pressure);
    }

    // Generate Brownian Kicks
    _brown->generateParallelKicks();

    // Set the starting guess for Lambda and K for
    // Augmented Lagrangian
    _constraint->setLagrangeCoeff(10.0);
    _constraint->setPenaltyCoeff(1000.0);

    // *************** Augmented Lagrangian Loop ************** //
    bool constraintMet = false;
    size_t alIter = 0, alMaxIter = 10;

    while (!constraintMet && (alIter < alMaxIter))
    {
      // Solve the unconstrained minimization
      _solver->solve();

      // Uzawa update
      _constraint->uzawaUpdate();

      // Update termination check quantities
      alIter++;
      constraintMet = _constraint->constraintSatisfied();
    }
    // *********************************************************//

    // Apply Kabsch Algorithm
    _ops->applyKabschAlgorithm();

    // Update kdTree, polyData and neighbors
    _ops->updatePolyData();
    _ops->updateNeighbors();

    // Calculate output values
    std::vector<double_t> msds(2, 0);
    msds = _ops->getMSD();

    // Write output to data file
    double_t morseEn = _ops->getMorseEnergy();
    double_t normEn = _ops->getNormalityEnergy();
    double_t circEn = _ops->getCircularityEnergy();
    double_t rmsAd = _ops->getRMSAngleDeficit();
    double_t volume = _ops->getVolume();
    double_t totalEn = morseEn + normEn + circEn;
    _detailedOP << _alpha << "\t" << _beta << "\t" << _gamma << "\t"
                << _ops->getAsphericity()
                << "\t"
                //<< morseEn << "\t"
                //<< normEn << "\t"
                //<< circEn << "\t"
                //<< _brown->getBrownianEnergy() << "\t"
                //<< _visco->getViscosityEnergy() << "\t"
                << totalEn << "\t" << _pressureBody->getPressureWork() << "\t"
                << msds[0] << "\t" << rmsAd << std::endl;
    // Store the energy and angular deficit data in circular buffer
    {
      // std::lock_guard<std::mutex> lock(_mut);
      _rmsAD->push_back(_step, rmsAd);
      _vol->push_back(_step, volume);
      _opsEn->push_back(_step, (morseEn + normEn + circEn));
    }

    // Update prevX
    _prevX = _x.head(3 * _N);

    // Send signal of step completion
    emit stepCompleted(_step++);
  }
}

vtkSmartPointer<vtkPolyData> LiveSimulation::GetPolyData()
{
  if (_computeVoronoi)
  {
    auto copyOpsPolyData = vtkSmartPointer<vtkPolyData>::New();
    copyOpsPolyData->DeepCopy(_ops->getPolyData());
    copyOpsPolyData->BuildLinks();
    size_t npts = copyOpsPolyData->GetNumberOfPoints();
    // get vertex positions
    std::vector<Vector3d> points(npts, Vector3d::Zero());
    // Calculate centroid of each triangle while updating points vector
    auto newPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = copyOpsPolyData->GetPolys();
    auto cellPointIds = vtkSmartPointer<vtkIdList>::New();
    cells->InitTraversal();
    while (cells->GetNextCell(cellPointIds))
    {
      size_t numCellPoints = cellPointIds->GetNumberOfIds();
      Vector3d centroid(0.0, 0.0, 0.0);
      for (size_t i = 0; i < numCellPoints; i++)
      {
        vtkIdType currCellPoint = cellPointIds->GetId(i);
        copyOpsPolyData->GetPoint(currCellPoint, &points[currCellPoint][0]);
        centroid += points[currCellPoint];
      }
      centroid /= numCellPoints;
      newPts->InsertNextPoint(&centroid[0]);
    }

    // Prepare valence Cell Data array
    auto valence = vtkSmartPointer<vtkIntArray>::New();
    valence->SetName("Valence");
    valence->SetNumberOfComponents(1);

    // Prepare new cell array for polygons
    auto newPolys = vtkSmartPointer<vtkCellArray>::New();
    for (size_t a = 0; a < npts; a++)
    {
      auto currPolyPtIds = vtkSmartPointer<vtkIdList>::New();
      std::list<neighbors> currPoly;
      Vector3d vec0, vecj, currCross, axis, centroid(0.0, 0.0, 0.0);
      double vec0_norm, vecj_norm, sign, currSin, currCos, currAngle;
      copyOpsPolyData->GetPointCells(a, currPolyPtIds);
      size_t numCellPoints = currPolyPtIds->GetNumberOfIds();
      // Get coordinates of first cell's centroid
      vtkIdType currId = currPolyPtIds->GetId(0);
      newPts->GetPoint(currId, &centroid[0]);
      vec0 = (centroid - points[a]).normalized();
      neighbors pt0(currId, 0.0);
      currPoly.push_back(pt0);
      // For remaining centroids
      for (auto j = 1; j < numCellPoints; j++)
      {
        currId = currPolyPtIds->GetId(j);
        newPts->GetPoint(currId, &centroid[0]);
        vecj = (centroid - points[a]).normalized();
        currSin = (vec0.cross(vecj)).norm();
        axis = (vec0.cross(vecj)).normalized();
        sign = axis.dot(points[a]);
        currSin = (sign > 0.0) ? currSin : -1.0 * currSin;
        currCos = vec0.dot(vecj);
        currAngle = (180 / M_PI) * atan2(currSin, currCos);
        currAngle = (currAngle < 0) ? (360 + currAngle) : currAngle;
        neighbors ptj(currId, currAngle);
        currPoly.push_back(ptj);
      }
      // Sort the list of neigbors and make a polygon
      currPoly.sort();
      newPolys->InsertNextCell(numCellPoints);
      for (auto t = currPoly.begin(); t != currPoly.end(); ++t)
      {
        neighbors n = *t;
        newPolys->InsertCellPoint(n._id);
      }
      valence->InsertNextTuple1(numCellPoints);
    }
    // Assign points and polygons to a new polydata and send it out
    auto newPolyData = vtkSmartPointer<vtkPolyData>::New();
    newPolyData->SetPoints(newPts);
    newPolyData->SetPolys(newPolys);
    newPolyData->GetCellData()->AddArray(valence);
    return newPolyData;
  }
  else
  {
    return _ops->getPolyData();
  }
}

void LiveSimulation::Reset()
{
  // Stop the simulation
  _keepRunning = false;
  // Set number of steps to 0
  _step = _resetStepVal;
  // Set positions and normals to starting value
  _x = _initialX;
  _prevX = _x.head(3 * _N);
  _ops->computeNormals();
  _ops->updatePolyData();
  _ops->updateNeighbors();
  // Clear the plot data
  _rmsAD->clear();
  _opsEn->clear();
  _vol->clear();
  // Reopen the output file discarding old content
  _detailedOP.close();
  _detailedOP.open(_dataOutputFile.c_str(), std::ofstream::out);
  _detailedOP << "Alpha"
              << "\t"
              << "Beta"
              << "\t"
              << "Gamma"
              << "\t"
              << "Asphericity"
              << "\t"
              //<< "MorseEn"
              //<< "\t"
              //<< "NormEn"
              //<< "\t"
              //<< "CircEn"
              //<< "\t"
              //<< "BrownEn"
              //<< "\t"
              //<< "ViscoEn"
              //<< "\t"
              << "TotalEn"
              << "\t"
              << "PressureWork"
              << "\t"
              << "MSD"
              << "\t"
              << "RMSAngleDeficit" << std::endl;
  emit updatePlotXAxis(_step);
  emit resetCompeleted();
}

void LiveSimulation::SaveScene(QString s)
{
  _ops->printVTKFile(s.toStdString());
}

} // namespace OPS
