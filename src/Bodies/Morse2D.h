#if !defined(__MORSE2D_H__)
#define __MORSE2D_H__

#include "Body.h"
#include "HelperFunctions.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cfloat>
#include <climits>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkVertexGlyphFilter.h>

namespace OPS {
// ************************* Morse2D Class ******************************//
/*!
 * \brief The Morse2D class
 * Computes energy and forces on a 2D lattice of particles with Morse
 */
class Morse2D : public Body {
public:
  typedef Eigen::Vector2d Vector2d;
  typedef Eigen::Matrix2Xd Matrix2Xd;
  typedef Eigen::Ref<Matrix2Xd> RefM2Xd;
  typedef Eigen::Map<Matrix2Xd> MapM2Xd;

  Morse2D(size_t n, double_t &f, RefM2Xd pos, RefM2Xd grad, double_t Lx,
          double_t Ly, double_t searchR);
  void applyKabschAlgorithm();
  virtual void compute();
  double_t getAverageEdgeLength();
  double_t getMeanSquaredDisplacement();
  double_t getMorseEnergy() { return _morseEn; }
  void printVTKFile(const std::string name);
  void reenterParticles();
  void saveInitialPosition() { _initialPositions = _positions; }
  void setMorseDistance(double_t r) { _re = r; }
  void setMorseWellWidth(double_t a) { _a = a; }
  void setSearchRadius(double_t s) { _searchRadius = s; }
  void setFVK(double_t g) { _gamma = g; }
  virtual void updateNeighbors();
  void updateDataForKabsch();

protected:
  double_t &_f, _msd = 0, _re = 1.0, _a = 4.62098120373, _gamma = 0.1,
                _searchRadius = 1.4, _morseEn;
  int _N;
  MapM2Xd _positions, _posGradient;
  Matrix2Xd _prevX, _initialPositions;
  std::vector<double_t> _bounds;
  std::vector<std::vector<int>> _boundaryCrossCounter;
  std::vector<std::vector<size_t>> _neighbors;
  std::vector<size_t> _initialNearestNeighbor;
};
//***************************************************************************//
} // namespace OPS
#endif //__MORSE2D_H__
