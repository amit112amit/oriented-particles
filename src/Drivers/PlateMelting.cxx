#include <stdio.h>
#include <string>
#include <vtkPolyDataReader.h>
#include "Morse2D.h"
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "ViscosityBody.h"

using namespace OPS;

int main(int argc, char* argv[]){
    clock_t t1, t2, t3;
    t1 = clock();

    if (argc != 2) {
	cout << "usage: " << argv[0] << " <filename>\n";
	return -1;
    }
    // ***************** Code flags *****************//
    bool loggingOn = false;

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];
    double_t Lx, Ly;

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->ReadAllVectorsOn();
    reader->Update();
    mesh = reader->GetOutput();
    double *bounds = mesh->GetBounds();
    Lx = bounds[1]-bounds[0];
    Ly = bounds[3]-bounds[2];
    // ********************************************************//

    // ******************* Read Simulation Parameters *********//

    double_t re=1.0, s=7.0, percentStrain = 15, searchR = 1.2, boxFactor = 1;
    size_t viterMax = 1000;
    size_t nameSuffix = 0;
    size_t step = 0;
    std::string baseFileName;

    InputParameters miscInp = OPS::readKeyValueInput( "miscInp.dat" );
    re = std::stod( miscInp["re"] );
    baseFileName = miscInp["baseFileName"];
    searchR = std::stod(miscInp["searchRadius"]);
    boxFactor = std::stod(miscInp["boxFactor"]);

    Lx += boxFactor*re;
    Ly += boxFactor*re;
    std::cout<< "Periodic box bounds : " << Lx << ", " << Ly << std::endl
	<< std::endl;

    // Input file should contain the following columns
    // Alpha Beta PercentStrain NumIterations PrintStep
    std::ifstream coolFile("schedule.dat");
    assert(coolFile);
    std::vector<std::vector<double_t> > coolVec;
    double_t currAlpha,currBeta,currPercentStrain,currViterMax,currPrintStep;

    std::string headerline;
    std::getline(coolFile, headerline);

    while (coolFile >> currAlpha >> currBeta >>  currPercentStrain
	    >> currViterMax >> currPrintStep) {
	std::vector<double> currLine;
	currLine.push_back(currAlpha);
	currLine.push_back(currBeta);
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

    // Read point coordinates from input mesh
    Eigen::Matrix2Xd coords(2,N);
    for(auto i = 0; i < N; ++i){
	double_t pt[3] = {0,0,0};
	mesh->GetPoint(i, pt);
	coords(0,i) = pt[0];
	coords(1,i) = pt[1];
    }

    // Prepare memory for energy, force
    double_t f;
    Eigen::VectorXd x(2*N), g(2*N), prevX(2*N);
    g.setZero(g.size());
    x.setZero(x.size());
    prevX.setZero(x.size());

    // Fill x with coords and rotVecs
    Eigen::Map<Eigen::Matrix2Xd> xpos(x.data(),2,N);
    xpos = coords;
    prevX = x.head(2*N);

    // Create Morse2D body
    Eigen::Map<Eigen::Matrix2Xd> posGrad(g.data(),2,N);
    Morse2D morseBody(N,f,xpos,posGrad,Lx,Ly,searchR);
    morseBody.setMorseDistance(re);
    s = 100*log(2.0)/(re*percentStrain);
    morseBody.setMorseWellWidth(s);

    // Create Brownian and Viscosity bodies
    Eigen::Map<Eigen::VectorXd> thermalX(x.data(),2*N,1);
    Eigen::Map<Eigen::VectorXd> thermalG(g.data(),2*N,1);
    double_t brownCoeff = 1.0, viscosity = 1.0;
    BrownianBody brown(2*N,brownCoeff,f,thermalX,thermalG,prevX);
    ViscosityBody visco(2*N,viscosity,f,thermalX,thermalG,prevX);

    // Create Model
    auto model = std::make_unique<Model>(2*N,f,g);
    model->addBody(std::make_shared<Morse2D>(morseBody));
    model->addBody(std::make_shared<BrownianBody>(brown));
    model->addBody(std::make_shared<ViscosityBody>(visco));
    // ****************************************************************//

    // ***************** Prepare Output Data files *********************//
    // Identify the Input structure name
    std::stringstream sstm;
    std::string dataOutputFile;
    std::string fname(baseFileName);

    // Detailed output data file
    ofstream detailedOP;
    sstm << fname << "-DetailedOutput.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    detailedOP.open(dataOutputFile.c_str(), std::ofstream::out);
    detailedOP << "#Alpha" << "\t"
	<< "Beta" << "\t"
	<< "MorseEn"  <<"\t"
	<< "BrownEn"  <<"\t"
	<< "ViscoEn"  <<"\t"
	<< "MSD"
	<< std::endl;

    // ************************* Create Solver ************************  //
    size_t m = 5, iprint = 1000, maxIter = 1e5;
    double_t factr = 10.0, pgtol = 1e-8;
    LBFGSBParams solverParams(m,iprint,maxIter,factr,pgtol);
    LBFGSBWrapper solver(solverParams, std::move(model), f, x, g);
    solver.turnOffLogging();
    // *****************************************************************//

    // ********************* Prepare data for simulation ****************//
    // Update the OPSBody member variables as per new positions
    morseBody.updateNeighbors();
    morseBody.saveInitialPosition(); /*!< For Mean Squared Displacement */
    // ******************************************************************//

    // ************************ OUTER SOLUTION LOOP **********************//
    size_t printStep;
    double_t alpha, beta;
    for(int z=0; z < coolVec.size(); z++){
	alpha = coolVec[z][0];
	beta = coolVec[z][1];
	percentStrain = coolVec[z][2];
	viterMax = coolVec[z][3];
	printStep = (int)coolVec[z][4];

	// Update OPS params
	s = (100 / (re*percentStrain))*log(2.0);
	morseBody.setMorseWellWidth(s);

	// For the very first iteration solve at zero temperature first
	if( z == 0 ){
	    brown.setCoefficient(0.0);
	    visco.setViscosity(0.0);
	    solver.solve();
	}

	// Update prevX
	prevX = x.head(2*N);

	// Set the viscosity and Brownian coefficient
	viscosity = alpha/(re*re);
	brownCoeff = std::sqrt( 2*alpha/beta )/re;
	if(loggingOn){
	    std::cout<< "Viscosity = " << viscosity << std::endl;
	    std::cout<< "Brownian Coefficient = " << brownCoeff
		<< std::endl;
	}
	brown.setCoefficient(brownCoeff);
	visco.setViscosity(viscosity);

	//**************  INNER SOLUTION LOOP ******************//
	for (int viter = 0; viter < viterMax; viter++) {
	    if(loggingOn)
		std::cout << std::endl
		    << "VISCOUS ITERATION: " << step
		    << std::endl
		    << std::endl;

	    //Generate Brownian kicks
	    brown.generateParallelKicks();

	    // Store data for Kabsch
	    //morseBody.updateDataForKabsch();

	    // Solve the unconstrained minimization
	    solver.solve();

	    // Apply Kabsch Algorithm
	    //morseBody.applyKabschAlgorithm();

	    // Re-enter the particles that left the box back into the box
	    morseBody.reenterParticles();

	    //Update neighbors
	    morseBody.updateNeighbors();

	    //********** Print relaxed configuration ************//
	    //We will print only after every currPrintStep iterations
	    if (viter % printStep == 0 && printStep <= viterMax) {
		sstm << fname << "-relaxed-" << nameSuffix++
		    <<".vtk";
		std::string rName = sstm.str();
		morseBody.printVTKFile(rName);
		sstm.str("");
		sstm.clear();
	    }

	    // Write output to data file
	    detailedOP
		<< alpha << "\t"
		<< beta << "\t"
		<< morseBody.getMorseEnergy() << "\t"
		<< brown.getBrownianEnergy() << "\t"
		<< visco.getViscosityEnergy() << "\t"
		<< morseBody.getMeanSquaredDisplacement()
		<< std::endl;

	    // Update prevX
	    prevX = x.head(2*N);
	    step++;
	}
    }
    // **********************************************************************//

    detailedOP.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
	<< " seconds" << std::endl;
    return 1;
}
