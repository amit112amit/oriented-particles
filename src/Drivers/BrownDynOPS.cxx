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

using namespace OPS;

int main(int argc, char* argv[]){
    clock_t t1, t2, t3;
    t1 = clock();

    if (argc != 2) {
	cout << "usage: " << argv[0] << " <filename>\n";
	return -1;
    }

    //******************** Optional parameters ********************//
    bool loggingOn = false;
    bool spikePrintOn = false;
    bool printAvgShapeData = false;

    // ***************** Read Input VTK File *****************//
    std::string inputFileName = argv[1];

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkSmartPointer<vtkPolyData> mesh;

    reader->SetFileName(inputFileName.c_str());
    reader->ReadAllVectorsOn();
    reader->Update();
    mesh = reader->GetOutput();
    // ********************************************************//

    // ******************* Read Simulation Parameters *********//

    // This flag determines if its a new simulation or a continuation
    bool continueFlag = false;

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
    continueFlag = std::stoi( miscInp["continueFlag"] );
    baseFileName = miscInp["baseFileName"];
    step = std::stoi( miscInp["step"] );
    nameSuffix = std::stoi( miscInp["nameSuffix"] );

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
    //Calculate the number of bonds;
    size_t numBonds = (int)((12*5 + (N-12)*6)/2);

    // Read point coordinates from input mesh
    Eigen::Matrix3Xd coords(3,N);
    for(auto i = 0; i < N; ++i){
	Eigen::Vector3d cp = Eigen::Vector3d::Zero();
	mesh->GetPoint(i, &(cp(0)));
	coords.col(i) = cp;
    }

    // Generate initial rotation vectors either from input point normals or
    // from starting point coordinates depending on continueFlag
    Eigen::Matrix3Xd rotVecs(3,N);
    if( continueFlag ){
	// Read normals from input file
	Eigen::Matrix3Xd normals(3,N);
	vtkSmartPointer< vtkDoubleArray > normalsArr =
	    vtkDoubleArray::FastDownCast(mesh->GetPointData()->GetNormals());
	normalsArr->vtkAbstractArray::SetNumberOfComponents( 3 );
	for(auto i = 0; i < N; ++i){
	    Eigen::Vector3d cp = Eigen::Vector3d::Zero();
	    normalsArr->GetTuple(i, &(cp(0)));
	    normals.col(i) = cp;
	}
	OPSBody::initialRotationVector(normals, rotVecs);
    }
    else{
	// Generate rotation vectors from input point coordinates
	OPSBody::initialRotationVector(coords, rotVecs);
    }

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
    OPSMesh ops(N,f,xpos,xrot,posGrad,rotGrad);
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
    std::stringstream sstm;
    std::string dataOutputFile;

    // Detailed output data file
    ofstream detailedOP;
    sstm << fname << "-DetailedOutput.dat";
    dataOutputFile = sstm.str();
    sstm.str("");
    sstm.clear();
    detailedOP.open(dataOutputFile.c_str(), std::ofstream::out |
	    std::ofstream::app);
    if( !continueFlag ){
	detailedOP << "#Step" <<"\t"
	    << "Alpha" << "\t"
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
	    << "MSD"
	    << std::endl;
    }

    // Create the output file for average data
    ofstream outerLoopFile;
    if( printAvgShapeData ){
	sstm << fname << "-AverageOutput.dat";
	dataOutputFile = sstm.str();
	sstm.str("");
	sstm.clear();
	outerLoopFile.open(dataOutputFile.c_str(), std::ofstream::out |
		std::ofstream::app);
	if( !continueFlag ){
	    outerLoopFile << "#BigStep" <<"\t"
		<< "PercentStrain" <<"\t"
		<< "Alpha" << "\t"
		<< "Beta" << "\t"
		<< "Gamma" << "\t"
		<< "PercentStrain" << "\t"
		<< "Radius"  <<"\t"
		<< "Asphericity"
		<< std::endl;
	}
    }

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
    if( !continueFlag ){
	// Renormalize positions such that avgEdgeLen = 1.0
	for(auto i=0; i < N; ++i){
	    xpos.col(i) = xpos.col(i)/avgEdgeLen;
	}
    }

//TODO:
// We cannot calculate initial position afresh when continuing a existing
// simulation
// We also cannot recalculate the average edge lengths because the structure
// would have distorted

    // Update the OPSBody member variables as per new positions
    ops.updatePolyData();
    ops.updateNeighbors();
    ops.saveInitialPosition(); /*!< For Mean Squared Displacement */
    avgEdgeLen = ops.getAverageEdgeLength();
    if(loggingOn)
	std::cout << "After renormalizing, Avg Edge Length = "
	    << avgEdgeLen << std::endl;

    // Calculate and set the spontaneous curvature
    double_t C0 = avgEdgeLen/ops.getAverageRadius();
    ops.setSpontaneousCurvature( C0 );
    // ******************************************************************//

    t3 = clock();
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
	s = (100 / (avgEdgeLen*percentStrain))*log(2.0);
	ops.setFVK(gamma);
	ops.setMorseWellWidth(s);

	// Set up the constraint value as the zero temperature value
	constraint->setConstraint(constrainedVal);

	// For the very first iteration solve at zero temperature first
	if( z == 0 && !continueFlag){
	    brown.setCoefficient(0.0);
	    visco.setViscosity(0.0);
	    solver.solve();
	}

	// Update prevX
	prevX = x.head(3*N);

	// Set the viscosity and Brownian coefficient
	viscosity = alpha/(avgEdgeLen*avgEdgeLen);
	brownCoeff = std::sqrt( 2*alpha/beta )/avgEdgeLen;
	if(loggingOn){
	    std::cout<< "Viscosity = " << viscosity << std::endl;
	    std::cout<< "Brownian Coefficient = " << brownCoeff
		<< std::endl;
	}
	brown.setCoefficient(brownCoeff);
	visco.setViscosity(viscosity);

	//**************  INNER SOLUTION LOOP ******************//
	Eigen::Matrix3Xd averagePosition( 3, N );
	averagePosition = Eigen::Matrix3Xd::Zero(3,N);

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

	    // Store data for Kabsch
	    ops.updateDataForKabsch();

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

	    // Add current solution to average position data
	    averagePosition += xpos;

	    // Calculate statistics about average displacement and
	    // largest displacement
	    double_t avgDisplacement, maxDisplacement;
	    Eigen::Matrix3Xd posDiff = xpos - prevPos;
	    avgDisplacement = posDiff.colwise().norm().sum()/N;
	    maxDisplacement = posDiff.colwise().norm().maxCoeff();
	    Eigen::Matrix3Xd normDiff(3,numBonds);
	    ops.getDiffNormals(normDiff);
	    Eigen::Matrix3Xd normals(3,numBonds);
	    ops.getNormals(normals);

	    //********** Print relaxed configuration ************//
	    //We will print only after every currPrintStep iterations
	    if (viter % printStep == 0 && printStep <= viterMax) {
		sstm << fname << "-relaxed-" << nameSuffix++
		    <<".vtk";
		std::string rName = sstm.str();
		ops.printVTKFile(rName);
		sstm.str("");
		sstm.clear();
	    }
	    // Print VTK file if there is an abrupt change in energy
	    if( spikePrintOn ){
		if ( std::abs( avgTotalEnergy ) > 0 &&
			std::abs((f - avgTotalEnergy)/avgTotalEnergy) > 3){
		    sstm << fname << "-Spike-" << step <<".vtk";
		    std::string rName = sstm.str();
		    ops.printVTKFile(rName);
		    sstm.str("");
		    sstm.clear();
		}
	    }

	    detailedOP << step << "\t"
		<< alpha << "\t"
		<< beta << "\t"
		<< gamma << "\t"
		<< ops.getAsphericity() << "\t"
		<< ops.getAverageRadius() << "\t"
		<< ops.getVolume() << "\t"
		<< ops.getArea() << "\t"
		<< ops.getMorseEnergy() << "\t"
		<< ops.getNormalityEnergy() << "\t"
		<< ops.getCircularityEnergy() << "\t"
		<< brown.getBrownianEnergy() << "\t"
		<< visco.getViscosityEnergy() << "\t"
		<< ops.getMeanSquaredDisplacement()
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

	if( printAvgShapeData ){
	    // Calculate the average particle positions and avg radius
	    double avgShapeRad = 0.0, avgShapeAsph = 0.0;
	    auto avgPos = vtkSmartPointer<vtkPoints>::New();
	    auto avgPosData = vtkSmartPointer<vtkDoubleArray>::New();
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
    }
    // *****************************************************************************//

    detailedOP.close();
    outerLoopFile.close();
    t2 = clock();
    float diff((float)t2 - (float)t1);
    std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
	<< " seconds" << std::endl;
    delete constraint;
    return 1;
}
