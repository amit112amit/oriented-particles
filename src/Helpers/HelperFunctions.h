#if !defined(__HELPERFUNCTIONS_H__)
#define __HELPERFUNCTIONS_H__

#include <iostream>
#include <random>
#include <stdlib.h>
#include <string>
#include <map>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkIdFilter.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace OPS{
typedef std::map< std::string, std::string > InputParameters;
typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix3Xd M3X;
typedef std::vector<vtkIdType> IdList;
typedef std::mt19937 Engine;
typedef std::normal_distribution<double_t> NormD;

class SimulationState{
public:
    SimulationState():N(0),nameSuffix(0),step(0),x(1),
        prevX(1),initPos(3,1),rng(0.,1.){}
    SimulationState(size_t n, size_t ns, size_t s, double_t g, double_t b,
	    Vec xi, Vec pi, M3X ip, IdList l, Engine e, NormD r): N(n),
        nameSuffix(ns), step(s), x(6*n), prevX(3*n), initPos(3,n),
        engine(e), rng(r), gamma(g), beta(b){
        nn = l; x = xi; prevX = pi; initPos = ip;
    }
    size_t getN(){return N;}
    size_t getNameSuffix(){return nameSuffix;}
    size_t getStep(){return step;}
    double_t getGamma(){return gamma;}
    double_t getBeta(){return beta;}
    Vec getX(){return x;}
    Vec getPrevX(){return prevX;}
    M3X getInitPos(){return initPos;}
    IdList getNeighbors(){return nn;}
    std::mt19937 getRandomEngine(){return engine;}
    NormD getRandomGenerator(){return rng;}
    void writeToFile(std::string file);
    static SimulationState readFromFile(std::string file);
private:
    size_t N;
    size_t nameSuffix;
    size_t step;
    double_t gamma = 0.01, beta = 10.0;
    static const size_t numsPerLine = 9;
    IdList nn;
    Vec x, prevX;
    M3X initPos;
    std::mt19937 engine;
    NormD rng;
};

//! Rotates a point cloud to get rid of rigid body motions
Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
                                      Eigen::Ref<Eigen::Matrix3Xd> out);

//! Generates a mesh for a topologically spherical point cloud
void delaunay3DSurf(vtkSmartPointer<vtkPolyData> poly, std::string fileName);

//! Calculate average edge length from a spherical point cloud
double_t getPointCloudAvgEdgeLen(std::string f);

//! Reads a key-value input file and returns a map object
InputParameters readKeyValueInput( std::string fileName );

//! Rotates a point cloud to get rid of rigid body motions
Eigen::Affine2d find2DAffineTransform(Eigen::Ref<Eigen::Matrix2Xd> in,
                                      Eigen::Ref<Eigen::Matrix2Xd> out);

//! Read Interpolation Data from Input File and return vector of vectors
void ReadInterpolationData(std::string, std::vector<double_t> &x,
                           std::vector<double_t> &r,
                           std::vector<double_t> &y,
                           std::vector<double_t> &e);

//! Return interpolated value based on data
std::vector<double_t> GetInterpolatedValue(double_t x,
                                           std::vector<double_t> &a,
                                           std::vector<double_t> &y,
                                           std::vector<double_t> &z,
                                           std::vector<double_t> &w);

void TestFind2DAffineTransform();
void TestFind3DAffineTransform();

}
#endif // __HELPERFUNCTIONS_H__
