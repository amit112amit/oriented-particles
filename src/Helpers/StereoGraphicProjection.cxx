#include "HelperFunctions.h"

namespace OPS {

Eigen::Matrix3Xd stereographicProjection(const Eigen::Matrix3Xd &p) {
  // Copy coordinates
  auto N = p.cols();
  Eigen::Matrix3Xd points(3, N);
  points = p;

  // Project points to unit sphere
  points.colwise().normalize();

  // Reset the center of the sphere to origin by translating
  points = points.colwise() - points.rowwise().mean();

  // Rotate all points so that the point in 0th column is along z-axis
  Eigen::Vector3d c = points.col(0);
  double_t cos_t = c(2);
  double_t sin_t = std::sqrt(1 - cos_t * cos_t);
  Eigen::Vector3d axis;
  axis << c(1), -c(0), 0.;
  Eigen::Matrix3d rotMat, axis_cross, outer;
  axis_cross << 0., -axis(2), axis(1), axis(2), 0., -axis(0), -axis(1), axis(0),
      0.;

  outer.noalias() = axis * axis.transpose();

  rotMat = cos_t * Eigen::Matrix3d::Identity() + sin_t * axis_cross +
           (1 - cos_t) * outer;
  Eigen::Matrix3Xd rPts(3, N);
  rPts = rotMat * points; // The points on a sphere rotated

  // Calculate the stereographic projections
  Eigen::Vector3d p0;
  Eigen::Map<Eigen::Matrix3Xd> l0(&(rPts(0, 1)), 3, N - 1);
  Eigen::Matrix3Xd l(3, N - 1), proj(3, N - 1);
  p0 << 0, 0, -1;
  c = rPts.col(0);
  l = (l0.colwise() - c).colwise().normalized();
  for (auto j = 0; j < N - 1; ++j) {
    proj.col(j) = ((p0(2) - l0(2, j)) / l(2, j)) * l.col(j) + l0.col(j);
  }
  return proj;
}

} // namespace OPS
