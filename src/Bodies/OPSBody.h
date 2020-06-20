/* Reference:
 *       Szeliski, Richard and Tonnesen, D. (1992). Surface Modeling with
 *       Oriented Particle Systems.	Siggraph ’92, 26(2), 160.
 *       https://doi.org/10.1017/CBO9781107415324.004
 */
#if !defined(__OPSBODY_H__)
#define __OPSBODY_H__

#include "Body.h"
#include "HelperFunctions.h"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <math.h>
#include <numeric>
#include <random>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkMassProperties.h>
#include <vtkOctreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

namespace OPS
{
// ************************* OPSBody Class *********************************//
/*!
 * \brief The OPSBody class
 * Computes energy and forces on a Oriented Particle System
 */
class OPSBody : public Body
{
public:
  typedef Eigen::Vector3d Vector3d;
  typedef Eigen::Map<Vector3d> MapV3d;
  typedef Eigen::Ref<const Vector3d> RefCV3d;
  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::Matrix3d Matrix3d;
  typedef Eigen::Ref<Matrix3d> RefM3d;
  typedef Eigen::Matrix3Xd Matrix3Xd;
  typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
  typedef Eigen::Map<Matrix3Xd> MapM3Xd;
  typedef Eigen::Matrix<double_t, 3, 4> Matrix3x4d;
  typedef Eigen::Matrix<double_t, 4, 3> Matrix4x3d;
  typedef Eigen::Quaterniond Quaterniond;
  typedef Eigen::AngleAxisd AngleAxisd;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point;
  typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
  typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
  typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
  typedef Delaunay::Vertex_handle Vertex_handle;
  typedef Delaunay::Face_circulator Face_circulator;

  OPSBody(size_t n, double_t &f, double_t R0, RefM3Xd pos, RefM3Xd rot,
          RefM3Xd posGrad, RefM3Xd rotGrad, RefM3Xd pX);
  void applyKabschAlgorithm();
  virtual void compute();
  void computeNormals();
  void diffNormalRotVec();
  double_t determineSearchRadius();
  double_t getAsphericity();
  double_t getArea();
  double_t getAverageEdgeLength();
  double_t getAverageRadius();
  double_t getCircularityEnergy() { return _circEn; }
  void getInitialPositions(Matrix3Xd &v) { v = _initialPositions; }
  double_t getMeanSquaredDisplacement();
  std::vector<double_t> getMSD();
  std::vector<size_t> getInitialNeighbors() { return _initialNearestNeighbor; }
  double_t getMorseEnergy() { return _morseEn; }
  double_t getNormalityEnergy() { return _normalEn; }
  vtkSmartPointer<vtkPolyData> getPolyData() { return _polyData; }
  double_t getRMSAngleDeficit();
  double_t getTotalEnergy();
  double_t getVolume();
  static void initialRotationVector(RefM3Xd pos, RefM3Xd rotVec);
  void printVTKFile(const std::string name);
  void resetToInitialPositions();
  void saveInitialPosition() { _initialPositions = _positions; }
  void setInitialNeighbors(std::vector<size_t> &x)
  {
    _initialNearestNeighbor = x;
  }
  void setInitialPositions(M3X p) { _initialPositions = p; }
  void setMorseDistance(double_t r) { _re = r; }
  void setMorseWellWidth(double_t a) { _a = a; }
  void setSearchRadius(double_t s) { _searchRadius = s; }
  void setFVK(double_t g) { _gamma = g; }
  void sphericalDelaunay();
  Delaunay stereoDelaunay();
  virtual void updateNeighbors();
  void updatePolyData();

protected:
  bool _updateArea = true, _updateRadius = true, _updateVolume = true;
  double_t &_f, _area = 0, _radius = 0, _volume = 0, _morseEn, _normalEn,
                _circEn, _msd_tgt = 0, _re = 1.0, _a = 4.62098120373,
                _gamma = 0.1, _searchRadius = 1.4, _msd;
  // The original zero temperature radius is stored here
  double_t _R0 = 0.0;
  int _N;
  MapM3Xd _positions, _posGradient, _rotGradient, _rotationVectors, _prevX;
  Matrix3Xd _normals, _initialPositions;
  std::vector<Matrix3d, Eigen::aligned_allocator<Matrix3d>> _diffNormalRV;
  std::vector<vtkSmartPointer<vtkIdList>> _neighbors;
  std::vector<size_t> _initialNearestNeighbor;
  vtkSmartPointer<vtkPolyData> _polyData;
  vtkSmartPointer<vtkOctreePointLocator> _octree;
  void updateRotationVectors();
};
//***************************************************************************//
} // namespace OPS
#endif //__OPSBODY_H__
