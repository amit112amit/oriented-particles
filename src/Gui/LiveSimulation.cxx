#include "LiveSimulation.h"

namespace OPS{

void LiveSimulation::Initialize(){
    // Initialize circular buffer with capacity 50000
    _rmsAD = new QwtCircBuffSeriesData(5000);
    _opsEn = new QwtCircBuffSeriesData(5000);
    _vol = new QwtCircBuffSeriesData(5000);

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = "T7.vtk";

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->ReadAllVectorsOn();
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//

    // ******************* Simulation Parameters *********//
    double_t re=1.0;
    double_t percentStrain = 15;
    double_t s = (100 / (re*percentStrain))*log(2.0);
    double_t brownCoeff = 1.0, viscosity = 1.0;
    ReadFvKAreaData("T7_OPS_Asphericity.dat");
    _zeroOpsEnVal = GetInterpolatedValue(_gamma,_opsEnDat);
    _zeroRmsAdVal = GetInterpolatedValue(_gamma,_rmsAdDat);
    _zeroVolVal = GetInterpolatedValue(_gamma,_volDat);

    // **********************************************************//

    // ***************** Create Bodies and Model ****************//
    // Set number of OPS particles
    _N = mesh->GetNumberOfPoints();

    // Read point coordinates from input mesh
    Eigen::Matrix3Xd coords(3,_N);
    for(auto i = 0; i < _N; ++i){
        Eigen::Vector3d cp = Eigen::Vector3d::Zero();
        mesh->GetPoint(i, &(cp(0)));
        coords.col(i) = cp;
    }

    // Generate initial rotation vectors either from starting point coordinates
    Eigen::Matrix3Xd rotVecs(3,_N);
    // Generate rotation vectors from input point coordinates
    OPSBody::initialRotationVector(coords, rotVecs);

    // Prepare memory for energy, force
    _x = Eigen::VectorXd(6*_N);
    _g = Eigen::VectorXd(6*_N);
    _prevX = Eigen::VectorXd(3*_N);
    _initialX = Eigen::VectorXd(6*_N);

    // Fill _x with coords and rotVecs
    Eigen::Map<Eigen::Matrix3Xd> xpos(_x.data(),3,_N), xrot(&(_x(3*_N)),3,_N),
            prevPos(_prevX.data(),3,_N);
    xpos = coords;
    xrot = rotVecs;
    _prevX = _x.head(3*_N);

    // Create OPSBody
    Eigen::Map<Eigen::Matrix3Xd> posGrad(_g.data(),3,_N), rotGrad(&_g(3*_N),3,_N);
    _ops = new OPSMesh(_N,_f,xpos,xrot,posGrad,rotGrad,prevPos);
    _ops->setFVK(_gamma);
    _ops->setMorseDistance(re);
    _ops->setMorseWellWidth(s);

    // Create InternalPressure body
    _pressureBody = new InternalPressure(_N,_f,xpos,prevPos,posGrad,_ops->getPolyData());

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(_x.data(),3*_N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(_g.data(),3*_N,1);
    _brown = new BrownianBody(3*_N,brownCoeff,_f,thermalX,thermalG,_prevX);
    _visco = new ViscosityBody(3*_N,viscosity,_f,thermalX,thermalG,_prevX);

    // Create area constraint
    vtkSmartPointer<vtkPolyData> poly = _ops->getPolyData();
    _constraint = new ExactAreaConstraint(_N, _f, xpos, posGrad, poly);
    _constraint->setConstraint( GetInterpolatedValue(_gamma,_areaDat) );

    // Create Model
    _model = new Model(6*_N,_f,_g);
    _model->addBody(_ops);
    _model->addBody(_brown);
    _model->addBody(_visco);
    _model->addBody(_constraint);
    _model->addBody(_pressureBody);
    // ****************************************************************//

    // ***************** Prepare Output Data files *********************//
    // Identify the Input structure name
    std::string fname = "T7";
    std::stringstream sstm;

    // Detailed output data file
    sstm << fname << "-DetailedOutput.dat";
    _dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    _detailedOP.open(_dataOutputFile.c_str(), std::ofstream::out);
    _detailedOP
                    << "Alpha" << "\t"
                    << "Beta" << "\t"
                    << "Gamma" << "\t"
                    << "Asphericity" << "\t"
                    //<< "MorseEn"  <<"\t"
                    //<< "NormEn"  <<"\t"
                    //<< "CircEn"  <<"\t"
                    //<< "BrownEn"  <<"\t"
                    //<< "ViscoEn"  <<"\t"
		    << "TotalEnergy" <<"\t"
                    << "PressureWork" <<"\t"
                    << "MSD" << "\t"
                    << "RMSAngleDeficit"
                    << std::endl;

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e5;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    _solver = new LBFGSBWrapper(solverParams, *(_model), _f, _x, _g);
    _solver->turnOffLogging();
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Calculate Average Edge Length
    double_t avgEdgeLen = _ops->getAverageEdgeLength();

    // Renormalize positions such that avgEdgeLen = 1.0
    for(auto i=0; i < _N; ++i){
        xpos.col(i) = xpos.col(i)/avgEdgeLen;
    }

    // Update the OPSBody member variables as per new positions
    _ops->updatePolyData();
    _ops->updateNeighbors();
    _ops->saveInitialPosition();
    avgEdgeLen = _ops->getAverageEdgeLength();

    // Relax at zero temperature once
    _brown->setCoefficient(0.0);
    _visco->setViscosity(0.0);
    _solver->solve();
    _prevX = _x.head(3*_N);

    // Save the current _x for resetting the simulation
    _initialX = _x;

    // Set Finite temperature coefficients
    viscosity = _alpha/(avgEdgeLen*avgEdgeLen);
    brownCoeff = std::sqrt( 2*_alpha/_beta )/avgEdgeLen;
    _brown->setCoefficient(brownCoeff);
    _visco->setViscosity(viscosity);
    // ******************************************************************//
    emit simulationReady();
}

// Interpolation of area
double_t LiveSimulation::GetInterpolatedValue(double_t x,
                                              std::vector<double_t> &y){
    int size = _gammaDat.size();
    int i = 0;
    if ( x >= _gammaDat[size - 2] ){
        i = size - 2;
    }
    else{
        while ( x > _gammaDat[i+1] ) i++;
    }
    double_t xL = _gammaDat[i], yL = y[i],
                    xR = _gammaDat[i+1], yR = y[i+1];
    if ( x < xL ) yR = yL;
    if ( x > xR ) yL = yR;
    double_t dydx = ( yR - yL ) / ( xR - xL );
    return yL + dydx * ( x - xL );
}

// Read gamma and area data from file
void LiveSimulation::ReadFvKAreaData(std::string s){
    ifstream inputFile(s.c_str());
    std::string line;
    double_t gamma, area, volume, rmsAd, opsEn, ignore;
    //Clean the existing vectors
    _gammaDat.clear();
    _volDat.clear();
    _areaDat.clear();
    _opsEnDat.clear();
    _rmsAdDat.clear();
    // Eat up the header line
    std::getline(inputFile, line);
    while (std::getline(inputFile, line)){
        std::istringstream ss(line);
        ss >> gamma >> ignore >> ignore >> volume
                        >> area >> ignore >> ignore >> ignore
                        >> opsEn >> rmsAd;
        _gammaDat.push_back(gamma);
        _volDat.push_back(volume);
        _areaDat.push_back(area);
        _opsEnDat.push_back(opsEn);
        _rmsAdDat.push_back(rmsAd);
    }
    //Reverse the vectors so that gamma values are increasing
    std::reverse(_gammaDat.begin(),_gammaDat.end());
    std::reverse(_volDat.begin(),_volDat.end());
    std::reverse(_areaDat.begin(),_areaDat.end());
    std::reverse(_opsEnDat.begin(),_opsEnDat.end());
    std::reverse(_rmsAdDat.begin(),_rmsAdDat.end());
}

// Update gamma and area constraint value
void LiveSimulation::UpdateGamma(double g){
    _gamma = g;
    _ops->setFVK(g);
    _constraint->setConstraint(GetInterpolatedValue(g,_areaDat));
    emit updateZeroOpsEn(GetInterpolatedValue(_gamma,_opsEnDat));
    emit updateZeroRmsAd(GetInterpolatedValue(_gamma,_rmsAdDat));
    emit updateZeroVolume(GetInterpolatedValue(_gamma,_volDat));
}

// Update pressure value
void LiveSimulation::UpdatePressure(double p){
    _pressure = p;
    _pressureBody->setPressure(p);
}

// Update beta and viscosity and brownian coefficients
void LiveSimulation::UpdateBeta(double b){
    if(b < 1e-6){
        _alpha = 0;
    }
    else{
        _alpha = 2.5e5;
        _beta = 1.0/b;
    }
    _brown->setCoefficient(std::sqrt(2*_alpha*b));
    _visco->setViscosity(_alpha);
}

// Start running
void LiveSimulation::SolveOneStep(){
    if(_keepRunning){
        if(_step == 0){
            // Relax at zero temperature once
            _brown->setCoefficient(0.0);
            _visco->setViscosity(0.0);
            _pressureBody->setPressure(0.0);
            _solver->solve();
            _prevX = _x.head(3*_N);
            _ops->updatePolyData();
            _ops->updateNeighbors();
            _brown->setCoefficient(std::sqrt(2*_alpha/_beta));
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

        while( !constraintMet && (alIter < alMaxIter)){
            // Solve the unconstrained minimization
            _solver->solve();

            //Uzawa update
            _constraint->uzawaUpdate();

            // Update termination check quantities
            alIter++;
            constraintMet = _constraint->constraintSatisfied();
        }
        // *********************************************************//

        // Apply Kabsch Algorithm
        _ops->applyKabschAlgorithm();

        //Update kdTree, polyData and neighbors
        _ops->updatePolyData();
        _ops->updateNeighbors();

        // Calculate output values
        std::vector<double_t> msds(2,0);
        msds = _ops->getMSD();

        // Write output to data file
        double_t morseEn = _ops->getMorseEnergy();
        double_t normEn = _ops->getNormalityEnergy();
        double_t circEn = _ops->getCircularityEnergy();
        double_t rmsAd = _ops->getRMSAngleDeficit();
        double_t volume = _ops->getVolume();
        double_t totalEn = morseEn + normEn + circEn;
        _detailedOP
                        << _alpha << "\t"
                        << _beta << "\t"
                        << _gamma << "\t"
                        << _ops->getAsphericity() << "\t"
                        //<< morseEn << "\t"
                        //<< normEn << "\t"
                        //<< circEn << "\t"
                        //<< _brown->getBrownianEnergy() << "\t"
                        //<< _visco->getViscosityEnergy() << "\t"
                        << totalEn <<"\t"
                        << _pressureBody->getPressureWork() <<"\t"
                        << msds[0] << "\t"
                        << rmsAd
                        << std::endl;
        // Store the energy and angular deficit data in circular buffer
        {
            //std::lock_guard<std::mutex> lock(_mut);
            _rmsAD->push_back(_step,rmsAd);
            _vol->push_back(_step,volume);
            _opsEn->push_back(_step,(morseEn+normEn+circEn));
        }

        // Update prevX
        _prevX = _x.head(3*_N);

        // Send signal of step completion
        emit stepCompleted(_step++);
    }
}

vtkSmartPointer<vtkPolyData> LiveSimulation::GetPolyData(){
    return _ops->getPolyData();
}

void LiveSimulation::Reset(){
    // Stop the simulation
    _keepRunning = false;
    // Set number of steps to 0
    _step = 0;
    // Set positions and normals to starting value
    _x = _initialX;
    _prevX = _x.head(3*_N);
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
    _detailedOP
                    << "Alpha" << "\t"
                    << "Beta" << "\t"
                    << "Gamma" << "\t"
                    << "Asphericity" << "\t"
                    << "MorseEn"  <<"\t"
                    << "NormEn"  <<"\t"
                    << "CircEn"  <<"\t"
                    << "BrownEn"  <<"\t"
                    << "ViscoEn"  <<"\t"
                    << "PressureWork" <<"\t"
                    << "MSD" << "\t"
                    << "RMSAngleDeficit"
                    << std::endl;
    emit resetCompeleted();
}

}
