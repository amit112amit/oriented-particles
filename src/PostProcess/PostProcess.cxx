#include <ctime>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <Eigen/Dense>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <boost/filesystem.hpp>
#include "SHTools.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned,K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;
typedef Delaunay::Face_circulator Face_circulator;
typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Point Point;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix3Xd Matrix3Xd;
typedef Eigen::Map<Matrix3Xd> Map3Xd;

int main(int argc, char* argv[]){

    // We expect three command line arguments:
    // 1. Number of particles in the simulation
    // 2. Path to the data file
    // 3. Output directory where output files will be created
    if( argc < 4 ){
	std::cout<< "Usage: ./postProcess <number_of_particles> <data_file> "
	    << "<output_dir>" << std::endl;
	return 0;
    }

    clock_t t = clock();

    // The number of particles in the simulation
    size_t N = std::stoi( std::string(argv[1]) );

    // Open the data file
    std::ifstream inputfile(argv[2]);
    if( !inputfile.is_open() ){
	std::cout<< "Unable to open file " << argv[2] << std::endl;
	return 0;
    }

    //**********************************************************************//
    // Create directory for storing the output files f.dat, Alm.dat and
    // spectrum.dat
    //

    auto outputDir = boost::filesystem::path( std::string( argv[3] ) );
    if( !boost::filesystem::create_directory( outputDir ) ){
	std::cout<< "Could not create output directory " << outputDir
	    << std::endl;
	return 0;
    }

    auto f_path = std::string( argv[3] ) + "/f.dat";
    std::ofstream f_file(f_path);
    auto Alm_path = std::string( argv[3] ) + "/Alm.dat";
    std::ofstream Alm_file(Alm_path);
    auto spec_path = std::string( argv[3] ) + "/spectrum.dat";
    std::ofstream spec_file(spec_path);

    if(!f_file){
	std::cout<< "Failed to create output file " << f_path;
	return 0;
    }
    if(!Alm_file){
	std::cout<< "Failed to create output file " << Alm_path;
	return 0;
    }
    if(!spec_file){
	std::cout<< "Failed to create output file " << spec_path;
	return 0;
    }

    //**********************************************************************//
    // Create the Gauss-Lengendre grid
    //
    int lmax = std::floor( std::sqrt(N) - 1 );
    int nlat, nlong;
    Eigen::VectorXd latglq(lmax + 1);
    Eigen::VectorXd longlq(2*lmax + 1);
    glqgridcoord_wrapper_(latglq.data(), longlq.data(), &lmax, &nlat, &nlong);
    auto numQ = nlat*nlong;

    Eigen::VectorXd gridglq((lmax + 1)*(2*lmax + 1));
    Eigen::VectorXd plx( (lmax + 1)*(lmax + 1)*(lmax + 2)/2 );
    Eigen::VectorXd w( lmax + 1 ), zero( lmax + 1 );
    Eigen::VectorXd cilm( 2*(lmax + 1)*(lmax + 1) );
    Eigen::VectorXd pspectrum( lmax + 1 );

    // Pre-compute all the matrices needed for expansion
    shglq_wrapper_(&lmax, zero.data(), w.data(), plx.data() );

    //Now we need to first identify which triangles each of the quadrature
    //points belongs to. For this, we need the coordinates of the points
    Matrix3Xd Q0(3,numQ);// Original quadrature points

    size_t index = 0;
    for( auto ip = 0; ip < nlong; ++ip){
	auto phi = longlq(ip)*M_PI/180.0;
	auto sin_p = std::sin(phi);
	auto cos_p = std::cos(phi);
	for( auto it = 0; it < nlat; ++it){
	    auto theta = (90.0 - latglq(it))*M_PI/180.0;
	    auto sin_t = std::sin(theta);
	    Q0.col(index++) << sin_t*cos_p, sin_t*sin_p, std::cos(theta);
	}
    }

    //**********************************************************************//
    // Some useful lambda functions
    //

    // A function to interpolate data to the quadrature points
    Matrix3Xd points(3,N);
    Matrix3Xd Q0_copy(3,numQ);
    VectorXd finput(N);
    auto interpolate = [&points, &Q0_copy, &gridglq, &finput](const size_t i,
	    const size_t j, const size_t k, const size_t q){
	Eigen::Vector3d v0 = points.col(i);
	Eigen::Vector3d v1 = points.col(j);
	Eigen::Vector3d v2 = points.col(k);
	Eigen::Vector3d qp = Q0_copy.col(q);
	auto A0 = ((v1-qp).cross((v2-qp))).norm();
	auto A1 = ((v2-qp).cross((v0-qp))).norm();
	auto A2 = ((v0-qp).cross((v1-qp))).norm();
	auto A = A0 + A1 + A2;
	gridglq(q) = (A0/A)*finput(i) + (A1/A)*finput(j) + (A2/A)*finput(k);
    };

    // Calculate the function to be expanded using Spherical Harmonics
    Eigen::VectorXd plms( (lmax + 1)*(lmax + 2)/2 );
    auto Plm = [&plms](const size_t l, const size_t m){
	return plms( l*(l+1)/2 + m );
    };

    // Calculate Alm from cilm
    auto A_lm = [&cilm, &lmax](const size_t l, const int m){
	int j = 0, n = m;
	if( m < 0 ){
	    j = 1;
	    n = -m;
	}
	return cilm( j + 2*(l + (lmax + 1)*n ) );
    };

    //**********************************************************************//
    // Now read the data file line by line
    //

    // The first line contains gamma and beta information
    std::string line, ignore;
    double_t gamma = 0, beta = 0;
    std::getline( inputfile, line );
    std::istringstream header( line );
    header >> ignore >> gamma >> ignore >> beta;

    // Read one row of data at a time. The data has positions and rotation
    // vectors. For the time being, we will ignore the rotation vectors.
    auto rowCount = 0;
    while( std::getline(inputfile, line) ){

	//std::cout<< "Processing row " << rowCount << std::endl;
	std::istringstream rowStream(line);

	// We will extract 3*N doubles representing particle positions
	// from the stream
	std::string value;
	size_t valCount = 0;
	Matrix3Xd origPoints(3,N);
	while( valCount < 3*N && rowStream.good() ){
	    std::getline(rowStream, value, ',');
	    origPoints( valCount%3, valCount/3 ) = std::stod(value);
	    valCount++;
	}

	// Calculate average radius and project points to average sphere
	auto R = origPoints.colwise().norm().mean();
	points = R*origPoints.colwise().normalized();
	Q0_copy = R*Q0.colwise().normalized();

	// Reset the center of the sphere to origin by translating
	Vector3d center = points.rowwise().mean();
	points = points.colwise() - center;

	// Calculate finput as a fucntion of theta and phi
	finput = points.colwise().norm() - R*Eigen::VectorXd::Ones(N);

	// Rotate all points so that the 0th point is along z-axis
	Vector3d c = points.col(0);
	double_t cos_t = c(2)/c.norm();
	double_t sin_t = std::sin( std::acos( cos_t ) );
	Vector3d axis;
	axis << c(1), -c(0), 0.;
	axis.normalize();
	Matrix3d rotMat, axis_cross, outer;
	axis_cross << 0. , -axis(2), axis(1),
		   axis(2), 0., -axis(0),
		   -axis(1), axis(0), 0.;

	outer.noalias() = axis*axis.transpose();

	rotMat = cos_t*Matrix3d::Identity() + sin_t*axis_cross +
	    (1 - cos_t)*outer;
	Matrix3Xd rPts(3,N);
	rPts = rotMat*points; // The points on a sphere rotated

	// Calculate the stereographic projections
	Vector3d p0;
	Map3Xd l0( &(rPts(0,1)), 3, N-1 );
	Matrix3Xd l(3,N-1), proj(3,N-1);
	p0 << 0,0,-R; // Point on the plane of projection
	c = rPts.col(0); // The point from which we are projecting
	l = (l0.colwise() - c).colwise().normalized(); // dirns of projections
	for( auto j=0; j < N-1; ++j ){
	    proj.col(j) = ((p0(2) - l0(2,j))/l(2,j))*l.col(j) + l0.col(j);
	}

	// Insert the projected points in a CGAL vertex_with_info vector
	std::vector< std::pair< Point, unsigned> > verts;
	for( auto j=0; j < N-1; ++j ){
	    verts.push_back(std::make_pair(Point(proj(0,j),proj(1,j)),j+1));
	}

	// Triangulate
	Delaunay dt( verts.begin(), verts.end() );

	// Rotate and project the quadrature points
	Matrix3Xd Qr(3,numQ);
	Qr = rotMat*Q0_copy; // Rotate the quadrature points

	Eigen::Matrix3Xd lQ(3,numQ), Qsp(3,numQ);
	lQ = (Qr.colwise() - c).colwise().normalized();//projn dirn unit vectors
	for( auto j=0; j < numQ; ++j ){
	    Qsp.col(j) = ((p0(2) - Qr(2,j))/lQ(2,j))*lQ.col(j) + Qr.col(j);
	}

	// Locate the quadrature points using stereographic triangulation
	for( auto j=0; j < numQ; ++j ){
	    auto query = Point( Qsp(0,j), Qsp(1,j) );
	    Delaunay::Locate_type lt;
	    int li;
	    Face_handle face = dt.locate( query, lt, li );
	    switch(lt){
		case Delaunay::FACE:
		    {
			auto id0 = face->vertex(0)->info();
			auto id1 = face->vertex(1)->info();
			auto id2 = face->vertex(2)->info();
			interpolate( id0, id1, id2, j );
			break;
		    }
		case Delaunay::EDGE:
		    {
			auto id1 = face->vertex( (li + 1)%3 )->info();
			auto id2 = face->vertex( (li + 2)%3 )->info();
			Eigen::Vector3d v1, v2, qp;
			v1 = points.col(id1);
			v2 = points.col(id2);
			qp = Q0_copy.col(j);
			double_t ratio = (qp - v1).norm()/(qp - v2).norm();
			gridglq(j) = (finput(id1) +
				ratio*finput(id2))/(1 + ratio);
			break;
		    }
		case Delaunay::VERTEX:
		    gridglq(j) = finput( face->vertex( li )->info() );
		    break;
		case Delaunay::OUTSIDE_CONVEX_HULL:
		    {
			Eigen::Vector3d v0, v1, v2;
			auto id0 = dt.is_infinite(face->vertex(0))?
			    0 : face->vertex(0)->info();
			auto id1 = dt.is_infinite(face->vertex(1))?
			    0 : face->vertex(1)->info();
			auto id2 = dt.is_infinite(face->vertex(2))?
			    0 : face->vertex(2)->info();
			interpolate( id0, id1, id2, j );
			break;
		    }
		default:
		    std::cout<< "Quadrature point " << j << " not found!"
			<< std::endl;
	    }
	}

	// Expand using Gauss-Legendre Quadrature
	shexpandglq_wrapper_( cilm.data(), &lmax, gridglq.data(), w.data(),
		plx.data());

	// Calculate Alm coefficients and their mean and variance
	index = 0;
	VectorXd Alm_vec( (lmax + 1)*(lmax + 1) );
	for( auto l = 0; l <= lmax; ++l ){
	    for( auto m = -l; m <= l; ++m){
		Alm_vec(index++) = A_lm(l,m);
	    }
	}
	Alm_file << Alm_vec.mean() << ","
	    << (Alm_vec - Alm_vec.mean()*
		    VectorXd::Ones((lmax + 1)*(lmax + 1)))
	    .array().square().mean() << std::endl;

	// Get the power spectrum and write it to file
	shpowerspectrum_wrapper_( cilm.data(), &lmax, pspectrum.data() );
	for(auto l = 0; l < lmax; ++l)
	    spec_file<< pspectrum(l) << ",";
	spec_file<< pspectrum(lmax + 1) << std::endl;

	// Calculate mean and variance of finput
	f_file << finput.mean() << ","
	    << (finput - finput.mean()*
		    VectorXd::Ones(N))
	    .array().square().mean() << std::endl;

	rowCount++;

    }
    inputfile.close();
    spec_file.close();
    f_file.close();
    Alm_file.close();

    std::cout<< "Time elapsed = " << ((float)(clock() - t))/CLOCKS_PER_SEC
	<< " seconds." << std::endl;

    return 1;
}
