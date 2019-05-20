#include "Morse2D.h"

namespace OPS {
// Constructor
Morse2D::Morse2D(size_t N, double_t &f, RefM2Xd pos, RefM2Xd grad, double_t Lx,
                 double_t Ly, double_t s)
    : _N(N), _f(f), _positions(pos.data(), 2, N),
      _posGradient(grad.data(), 2, N), _searchRadius(s) {
  // Initialize internal arrays
  _prevX = _positions;

  // Initialize the vectors
  _bounds.push_back(Lx);
  _bounds.push_back(Ly);
  _boundaryCrossCounter =
      std::vector<std::vector<int>>(_N, std::vector<int>(2, 0));
  _neighbors = std::vector<std::vector<size_t>>(_N, std::vector<size_t>(0));
  // Calculate neighbors
  for (auto i = 0; i < _N; ++i) {
    std::vector<size_t> currNeigh;
    for (auto j = 0; j < _N; ++j) {
      if (j != i) {
        Vector2d dX = _positions.col(i) - _positions.col(j);
        for (int k = 0; k < 2; ++k) {
          if (dX(k) > _bounds[k] * 0.5)
            dX(k) -= _bounds[k];
          if (dX(k) <= -1 * _bounds[k] * 0.5)
            dX(k) += _bounds[k];
        }
        if (dX.norm() < _searchRadius) {
          currNeigh.push_back(j);
        }
      }
    }
    _neighbors[i] = currNeigh;
  }
  // Set the initial nearest neighbor vector
  for (auto i = 0; i < _N; ++i) {
    double_t nearestDist = DBL_MAX;
    size_t nearest = UINT_MAX;
    for (auto j : _neighbors[i]) {
      double_t dist = (_positions.col(i) - _positions.col(j)).norm();
      if (dist < nearestDist) {
        nearestDist = dist;
        nearest = j;
      }
    }
    _initialNearestNeighbor.push_back(nearest);
  }
  updateNeighbors();
  return;
}

// Move particles outside the box, back in the box. We assume that the box
// is centered at the origin and is symmetric about the x and y axes
void Morse2D::reenterParticles() {
  for (auto i = 0; i < _N; ++i) {
    Vector2d pos = _positions.col(i);
    for (auto j = 0; j < 2; ++j) {
      if (pos(j) <= -0.5 * _bounds[j]) {
        std::cout << "Particle " << i << " crossed " << j << " low "
                  << std::endl;
        _positions(j, i) += _bounds[j];
        _boundaryCrossCounter[i][j]--;
      }
      if (pos(j) > 0.5 * _bounds[j]) {
        std::cout << "Particle " << i << " crossed " << j << " high "
                  << std::endl;
        _positions(j, i) -= _bounds[j];
        _boundaryCrossCounter[i][j]++;
      }
    }
  }
}

// Update neighbors
void Morse2D::updateNeighbors() {
  for (auto i = 0; i < _N; ++i) {
    std::vector<size_t> currNeigh;
    for (auto j = 0; j < _N; ++j) {
      if (j != i) {
        Vector2d dX = _positions.col(i) - _positions.col(j);
        for (int k = 0; k < 2; ++k) {
          if (dX(k) > _bounds[k] * 0.5)
            dX(k) -= _bounds[k];
          if (dX(k) <= -1 * _bounds[k] * 0.5)
            dX(k) += _bounds[k];
        }
        if (dX.norm() < _searchRadius) {
          currNeigh.push_back(j);
        }
      }
    }
    _neighbors[i] = currNeigh;
  }
  // Remove i as a neighbor from list of j's neighbors to avoid
  // double counting
  for (auto i = 0; i < _N; ++i) {
    for (auto j : _neighbors[i]) {
      _neighbors[j].erase(
          std::remove(_neighbors[j].begin(), _neighbors[j].end(), i),
          _neighbors[j].end());
    }
  }
}

// Print VTK file
void Morse2D::printVTKFile(const std::string fileName) {
  vtkNew<vtkPoints> pts;
  vtkNew<vtkPolyData> poly, out;
  vtkNew<vtkVertexGlyphFilter> vgf;
  vtkNew<vtkPolyDataWriter> writer;
  for (auto i = 0; i < _N; ++i) {
    pts->InsertNextPoint(_positions(0, i), _positions(1, i), 0);
  }
  poly->SetPoints(pts.Get());
  vgf->SetInputData(poly.Get());
  vgf->Update();
  writer->SetInputConnection(vgf->GetOutputPort());
  writer->SetFileName(fileName.c_str());
  writer->Write();
}

// Update data for Kabsch algorithm
void Morse2D::updateDataForKabsch() { _prevX = _positions; }

// Apply Kabsch algorithm
void Morse2D::applyKabschAlgorithm() {
  Eigen::Affine2d A;
  A = find2DAffineTransform(_positions, _prevX);
  if (A.translation().norm() > 4)
    std::cout << "Kabsch translation = " << A.translation().transpose()
              << std::endl;
  for (auto i = 0; i < _N; ++i) {
    _positions.col(i) = A.linear() * _positions.col(i) + A.translation();
  }
}

// Calculate Mean Squared displacement
double_t Morse2D::getMeanSquaredDisplacement() {
  _msd = 0;
  for (auto i = 0; i < _N; ++i) {
    Vector2d xi, xj;
    xi = _positions.col(i) - _initialPositions.col(i);
    xj = _positions.col(_initialNearestNeighbor[i]) -
         _initialPositions.col(_initialNearestNeighbor[i]);
    for (auto k = 0; k < 2; ++k) {
      xi(k) += _boundaryCrossCounter[i][k] * _bounds[k];
      xj(k) +=
          _boundaryCrossCounter[_initialNearestNeighbor[i]][k] * _bounds[k];
    }
    _msd += (xi - xj).squaredNorm();
  }
  // We have double counted the displacements. So we need to divide by 2
  _msd /= 2 * _N;
  return _msd;
}

// Compute function
void Morse2D::compute() {
  _morseEn = 0.0;
  for (auto i = 0; i < _N; ++i) {
    for (auto j : _neighbors[i]) {
      double_t morse, exp_1, exp_2, r;
      Vector2d dX, dMdr;
      dX = _positions.col(j) - _positions.col(i);
      for (auto k = 0; k < 2; ++k) {
        if (dX(k) > _bounds[k] * 0.5)
          dX(k) -= _bounds[k];
        if (dX(k) <= -1 * _bounds[k] * 0.5)
          dX(k) += _bounds[k];
      }
      r = dX.norm();
      exp_1 = exp(-_a * (r - _re));
      exp_2 = exp_1 * exp_1;
      morse = exp_2 - 2 * exp_1;
      dMdr = (2 * _a / r) * (exp_1 - exp_2) * dX;
      _morseEn += morse;
      _posGradient.col(i) -= dMdr;
      _posGradient.col(i) += dMdr;
    }
  }
  return;
}
} // namespace OPS
