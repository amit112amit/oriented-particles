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
typedef Eigen::Map<Matrix3Xd> MapM3Xd;
typedef Eigen::Ref<Matrix3Xd> RefM3Xd;

int main(int argc, char* argv[]){

    // We expect three command line arguments:
    // 1. Path to the data file
    // 2. Output directory where output files will be created
    if( argc < 3 ){
	std::cout<< "Usage: ./postProcess <data_file> <output_dir>"
	    << std::endl;
	return 0;
    }

    clock_t t = clock();

    // Open the data file
    std::ifstream inputfile(argv[1]);
    if( !inputfile.is_open() ){
	std::cout<< "Unable to open file " << argv[1] << std::endl;
	return 0;
    }

    //**********************************************************************//
    // Create directory for storing the output files f.dat, Alm.dat and
    // spectrum.dat
    //

    auto outputDir = boost::filesystem::path( std::string( argv[2] ) );
    if( !boost::filesystem::create_directory( outputDir ) ){
	std::cout<< "Could not create output directory " << outputDir
	    << std::endl;
	return 0;
    }

    auto f_path = std::string( argv[2] ) + "/f.dat";
    std::ofstream f_file(f_path);
    auto u_path = std::string( argv[2] ) + "/u2.dat";
    std::ofstream u_file(u_path);
    auto Alm_path = std::string( argv[2] ) + "/Alm.dat";
    std::ofstream Alm_file(Alm_path);
    auto spec_path = std::string( argv[2] ) + "/spectrum.dat";
    std::ofstream spec_file(spec_path);

    if(!f_file){
	std::cout<< "Failed to create output file " << f_path;
	return 0;
    }
    if(!u_file){
	std::cout<< "Failed to create output file " << u_path;
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
    // Now read the first( or the header ) line of the data file
    // The first line contains N, gamma and beta information
    std::string line, ignore;
    size_t N;
    double_t gamma = 0, beta = 0;
    std::getline( inputfile, line );
    std::istringstream header( line );
    header >> ignore >> N >> ignore >> gamma >> ignore >> beta;

    //**********************************************************************//
    // Create the Gauss-Lengendre grid
    //
    int lmax = std::floor( std::sqrt(N) - 1 );
    int nlat, nlong;
    Eigen::VectorXd latglq(lmax + 1);
    Eigen::VectorXd longlq(2*lmax + 1);
    glqgridcoord_wrapper_(latglq.data(), longlq.data(), &lmax, &nlat, &nlong);
    auto numQ = nlat*nlong;

    Eigen::VectorXd gridglq(numQ);
    Eigen::VectorXd plx( (lmax + 1)*(lmax + 1)*(lmax + 2)/2 );
    Eigen::VectorXd w( lmax + 1 ), zero( lmax + 1 );
    Eigen::VectorXd cilm( 2*(lmax + 1)*(lmax + 1) );
    Eigen::VectorXd pspectrum( lmax + 1 );

    // Pre-compute all the matrices needed for expansion
    shglq_wrapper_(&lmax, zero.data(), w.data(), plx.data() );

    //**********************************************************************//
    // Some useful lambda functions
    //

    // A function to read N rows of a 3xN matrix from a stringstream
    auto readMatrix3Xd = [](std::istringstream &ss, RefM3Xd X){
	auto N = X.cols();
	size_t valCount = 0;
	while( valCount < 3*N && ss.good() ){
	    std::string value;
	    std::getline(ss, value, ',');
	    X(valCount%3, valCount/3) = std::stod(value);
	    valCount++;
	}
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

    // The first row of data is ALWAYS the zero temperature particle position
    // and rotation vectors. We will use it to locate the triangles containing
    // the quadrature points for interpolation and also the zero temperature
    // radius
    Matrix3Xd zeroTempX(3,N);
    std::getline(inputfile, line);
    std::istringstream zeroTempRow(line);
    readMatrix3Xd( zeroTempRow, zeroTempX );

    // Calculate the zero temperature radius
    double_t R0 = zeroTempX.colwise().norm().mean();

    // Now we need to project the zeroTempX to a sphere of radius R0
    Matrix3Xd X0(3,N);
    X0 = R0*zeroTempX.colwise().normalized();

    // The following matrix is unit normals associated with X0
    Matrix3Xd x0(3,N);
    x0 = X0.normalized();

    // Quadrature point coordinates in Cartesian on a sphere of radius R0
    Matrix3Xd Q0(3,numQ);
    size_t index = 0;
    for( auto ip = 0; ip < nlong; ++ip){
	auto phi = longlq(ip)*M_PI/180.0;
	auto sin_p = std::sin(phi);
	auto cos_p = std::cos(phi);
	for( auto it = 0; it < nlat; ++it){
	    auto theta = (90.0 - latglq(it))*M_PI/180.0;
	    auto sin_t = std::sin(theta);
	    Q0.col(index++) << R0*sin_t*cos_p, R0*sin_t*sin_p,
		R0*std::cos(theta);
	}
    }

    // Rotate all points of the shell so that the 0th point is along z-axis
    Vector3d c = X0.col(0);
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
    rotMat = cos_t*Matrix3d::Identity() + sin_t*axis_cross + (1 - cos_t)*outer;
    Matrix3Xd rPts(3,N);
    rPts = rotMat*X0;

    // Calculate the stereographic projections
    Vector3d p0;
    p0 << 0,0,-R0; // Point on the plane of projection
    c = rPts.col(0); // The point from which we are projecting

    MapM3Xd l0( &(rPts(0,1)), 3, N-1 );
    Matrix3Xd l(3,N-1), proj(3,N-1);
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

    // Rotate the quadrature points
    Matrix3Xd Qr(3,numQ);
    Qr = rotMat*Q0;

    // Stereographic projection of quadrature points
    Eigen::Matrix3Xd lQ(3,numQ), Qsp(3,numQ);
    lQ = (Qr.colwise() - c).colwise().normalized();
    for( auto j=0; j < numQ; ++j ){
	Qsp.col(j) = ((p0(2) - Qr(2,j))/lQ(2,j))*lQ.col(j) + Qr.col(j);
    }
    //**********************************************************************//

    //******************** Determine interpolation weights ******************//
    std::vector< std::vector< std::pair<size_t, double_t> > > interpolData;
    for( auto q=0; q < numQ; ++q ){
	auto query = Point( Qsp(0,q), Qsp(1,q) );
	Delaunay::Locate_type lt;
	int li;
	Face_handle face = dt.locate( query, lt, li );
	std::vector< std::pair< size_t, double_t > > nodesAndWeights;
	switch(lt){
	    case Delaunay::FACE:
		{
		    auto i = face->vertex(0)->info();
		    auto j = face->vertex(1)->info();
		    auto k = face->vertex(2)->info();
		    Eigen::Vector3d v0 = X0.col(i);
		    Eigen::Vector3d v1 = X0.col(j);
		    Eigen::Vector3d v2 = X0.col(k);
		    Eigen::Vector3d qp = Q0.col(q);
		    auto A0 = ((v1-qp).cross((v2-qp))).norm();
		    auto A1 = ((v2-qp).cross((v0-qp))).norm();
		    auto A2 = ((v0-qp).cross((v1-qp))).norm();
		    auto A = A0 + A1 + A2;
		    nodesAndWeights.push_back(
			    std::make_pair(i, A0/A));
		    nodesAndWeights.push_back(
			    std::make_pair(j, A1/A));
		    nodesAndWeights.push_back(
			    std::make_pair(k, A2/A));
		    break;
		}
	    case Delaunay::OUTSIDE_CONVEX_HULL:
		{
		    auto i = dt.is_infinite(face->vertex(0))?
			0 : face->vertex(0)->info();
		    auto j = dt.is_infinite(face->vertex(1))?
			0 : face->vertex(1)->info();
		    auto k = dt.is_infinite(face->vertex(2))?
			0 : face->vertex(2)->info();
		    Eigen::Vector3d v0 = X0.col(i);
		    Eigen::Vector3d v1 = X0.col(j);
		    Eigen::Vector3d v2 = X0.col(k);
		    Eigen::Vector3d qp = Q0.col(q);
		    auto A0 = ((v1-qp).cross((v2-qp))).norm();
		    auto A1 = ((v2-qp).cross((v0-qp))).norm();
		    auto A2 = ((v0-qp).cross((v1-qp))).norm();
		    auto A = A0 + A1 + A2;
		    nodesAndWeights.push_back(
			    std::make_pair(i, A0/A));
		    nodesAndWeights.push_back(
			    std::make_pair(j, A1/A));
		    nodesAndWeights.push_back(
			    std::make_pair(k, A2/A));
		    break;
		}
	    case Delaunay::EDGE:
		{
		    auto id1 = face->vertex( (li + 1)%3 )->info();
		    auto id2 = face->vertex( (li + 2)%3 )->info();
		    Eigen::Vector3d v1, v2, qp;
		    v1 = X0.col(id1);
		    v2 = X0.col(id2);
		    qp = Q0.col(q);
		    double_t ratio = (qp - v1).norm()/(qp - v2).norm();
		    nodesAndWeights.push_back(
			    std::make_pair(id1, 1.0/(1 + ratio)));
		    nodesAndWeights.push_back(
			    std::make_pair(id2, ratio/(1 + ratio)));
		    break;
		}
	    case Delaunay::VERTEX:
		nodesAndWeights.push_back(std::make_pair(li, 1.0));
		break;
	    default:
		std::cout<< "Quadrature point " << q << " not found!"
		    << std::endl;
	}
	interpolData.push_back( nodesAndWeights );
    }
    //**********************************************************************//

    //***************************** MAIN LOOP ******************************//
    // Now iterate over the remaining rows
    auto rowCount = 1;
    while( std::getline(inputfile, line) ){

	std::istringstream rowStream(line);

	// We will extract 3*N doubles representing particle positions
	// from the stream and ignore the 3*N rotation vector components
	Matrix3Xd Xt(3,N);
	readMatrix3Xd( rowStream, Xt );

	// Calculate finput for all particles
	VectorXd finput(N);
	finput = ((Xt - X0).array()*x0.array()).matrix().colwise().sum();
	auto f0 = finput.mean();

	// Calculate tangent fluctuations
	VectorXd u_squared(N);
	for( auto i = 0; i < N; ++i){
	    Vector3d u, ft;
	    u = x0.col(i).cross( Xt.col(i).cross( x0.col(i) ) );
	    u_squared(i) = u.squaredNorm();
	}

	// Now calculate f at quadrature points from finput by interpolation
	gridglq.setZero(numQ);
	for(auto q = 0; q < numQ; ++q){
	    auto nodeWeights = interpolData[q];
	    for(const auto &nw : nodeWeights){
		gridglq(q) += finput(nw.first)*nw.second;
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
	spec_file<< pspectrum(lmax) << std::endl;

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
    u_file.close();
    Alm_file.close();

    std::cout<< "Time elapsed = " << ((float)(clock() - t))/CLOCKS_PER_SEC
	<< " seconds." << std::endl;

    return 1;
}
