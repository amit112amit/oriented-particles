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
	std::cout<< "Usage: ./calcVolume <data_file> <output_dir>"
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
    //

    auto outputFile = std::string( argv[2] ) + "/volume.dat";
    std::ofstream vol_file( outputFile );

    if(!vol_file){
	std::cout<< "Failed to create output file " << outputFile << std::endl;
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
    int lmax = std::floor( std::sqrt(N) - 1 );

    // Write column names in the output files
    vol_file << "volume" << std::endl;
    vol_file.precision(7);

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
    //**********************************************************************//

    //***************************** MAIN LOOP ******************************//
    // Skip the first row which is just zero temperature data
    std::getline( inputfile, line );

    // Now iterate over the remaining rows
    auto rowCount = 1;
    while( std::getline(inputfile, line) ){

	std::istringstream rowStream(line);

	// We will extract 3*N doubles representing particle positions
	// from the stream and ignore the 3*N rotation vector components
	Matrix3Xd Xt(3,N), X0(3,N);
	readMatrix3Xd( rowStream, Xt );

	// Project the points to a unit sphere
	X0 = Xt.colwise().normalized();

	// Rotate all points of the shell so that the 0th point is along z-axis
	Vector3d c = X0.col(0);
	double_t cos_t = c(2);
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
	p0 << 0,0,-1.0; // Point on the plane of projection
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

	// Iterate over the triangles to calculate volume
	auto volume = 0.0;
	for( auto ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end();
	       	++ffi){
	    auto i = ffi->vertex(0)->info();
	    auto j = ffi->vertex(2)->info();
	    auto k = ffi->vertex(1)->info();
	    volume += 0.166666667*(Xt.col(i).dot(Xt.col(j).cross(Xt.col(k))));
	}

	// Iterate over infinite faces
	Face_circulator fc = dt.incident_faces(dt.infinite_vertex()), done(fc);
	if (fc != 0) {
	    do{
		auto i = dt.is_infinite(fc->vertex(0))?0:fc->vertex(0)->info();
		auto j = dt.is_infinite(fc->vertex(2))?0:fc->vertex(2)->info();
		auto k = dt.is_infinite(fc->vertex(1))?0:fc->vertex(1)->info();
		volume += 0.166666667*
		    (Xt.col(i).dot(Xt.col(j).cross(Xt.col(k))));
	    }while(++fc != done);
	}
	vol_file << volume << std::endl;
	rowCount++;
    }
    inputfile.close();
    vol_file.close();

    std::cout<< "Time elapsed = " << ((float)(clock() - t))/CLOCKS_PER_SEC
	<< " seconds." << std::endl;

    return 1;
}
