#if !defined(__HELPERFUNCTIONS_H__)
#define __HELPERFUNCTIONS_H__

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <map>
#include <random>
#include <stdlib.h>
#include <string>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

namespace OPS
{
typedef std::map<std::string, std::string> InputParameters;
typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix3Xd M3X;
typedef Eigen::Ref<M3X> RefM3X;
typedef std::vector<size_t> IdList;
typedef std::mt19937 Engine;
typedef std::normal_distribution<double_t> NormD;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef Delaunay::Face_circulator Face_circulator;

class SimulationState
{
public:
  SimulationState()
      : N(0), nameSuffix(0), step(0), x(1), prevX(1), initPos(3, 1),
        rng(0., 1.), radius0(1.0) {}
  SimulationState(size_t n, size_t ns, size_t s, double_t g, double_t b, double_t rad,
                  Vec xi, Vec pi, M3X ip, IdList l, Engine e, NormD r)
      : N(n), nameSuffix(ns), step(s), x(6 * n), prevX(3 * n), initPos(3, n),
        engine(e), rng(r), gamma(g), beta(b), radius0(rad)
  {
    nn = l;
    x = xi;
    prevX = pi;
    initPos = ip;
  }
  size_t getN() { return N; }
  size_t getNameSuffix() { return nameSuffix; }
  size_t getStep() { return step; }
  double_t getGamma() { return gamma; }
  double_t getBeta() { return beta; }
  double_t getRadius0() { return radius0; }
  Vec getX() { return x; }
  Vec getPrevX() { return prevX; }
  M3X getInitPos() { return initPos; }
  IdList getNeighbors() { return nn; }
  std::mt19937 getRandomEngine() { return engine; }
  NormD getRandomGenerator() { return rng; }
  void writeToFile(std::string file);
  static SimulationState readFromFile(std::string file);

private:
  size_t N;
  size_t nameSuffix;
  size_t step;
  double_t gamma = 0.01, beta = 10.0, radius0 = 1.0;
  static const size_t numsPerLine = 9;
  IdList nn;
  Vec x, prevX;
  M3X initPos;
  std::mt19937 engine;
  NormD rng;
};

//! A struct to store a vtkIdType and an angle
struct neighbors
{
  vtkIdType _id;
  double_t _angle;
  neighbors(vtkIdType i, double_t a) : _id(i), _angle(a) {}
  bool operator<(const neighbors &n) const { return _angle < n._angle; }
};

//! Rotates a point cloud to get rid of rigid body motions
Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
                                      Eigen::Ref<Eigen::Matrix3Xd> out);

//! Generates a mesh for a topologically spherical point cloud
void delaunay3DSurf(vtkSmartPointer<vtkPolyData> poly, std::string fileName);

//! Generates a mesh for a topologically spherical point cloud
Delaunay delaunayStereo(Eigen::Ref<Eigen::Matrix3Xd> p);

//! Calculate average edge length from a spherical point cloud
double_t getPointCloudAvgEdgeLen(std::string f);

//! Reads a key-value input file and returns a map object
InputParameters readKeyValueInput(std::string fileName);

//! Rotates a point cloud to get rid of rigid body motions
Eigen::Affine2d find2DAffineTransform(Eigen::Ref<Eigen::Matrix2Xd> in,
                                      Eigen::Ref<Eigen::Matrix2Xd> out);

//! Read Interpolation Data from Input File and return vector of vectors
void ReadInterpolationData(std::string, std::vector<double_t> &x,
                           std::vector<double_t> &r, std::vector<double_t> &y,
                           std::vector<double_t> &e);

//! Return interpolated value based on data
std::vector<double_t> GetInterpolatedValue(double_t x, std::vector<double_t> &a,
                                           std::vector<double_t> &y,
                                           std::vector<double_t> &z,
                                           std::vector<double_t> &w);

//! Stereographic projection approach to mesh shells
Eigen::Matrix3Xd stereographicProjection(const Eigen::Matrix3Xd &p);

void TestFind2DAffineTransform();
void TestFind3DAffineTransform();

} // namespace OPS
#endif // __HELPERFUNCTIONS_H__
