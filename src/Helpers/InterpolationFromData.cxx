#include "HelperFunctions.h"

using namespace std;
typedef std::vector<double_t> vec;

namespace OPS {

// Read gamma and area data from file
void ReadInterpolationData(string s, vec &gammaDat, vec &volDat, vec &energyDat,
                           vec &radData) {
  ifstream inputFile(s.c_str());
  string line;
  double_t gamma, volume, energy, radius, ignore;
  // Clean the existing vectors
  gammaDat.clear();
  volDat.clear();
  energyDat.clear();
  radData.clear();
  // Eat up the header line
  getline(inputFile, line);
  while (getline(inputFile, line)) {
    istringstream ss(line);
    ss >> gamma >> ignore >> radius >> volume >> ignore >> ignore >> ignore >>
        ignore >> energy >> ignore;
    gammaDat.push_back(gamma);
    volDat.push_back(volume);
    energyDat.push_back(energy);
    radData.push_back(radius);
  }
  // Reverse the vectors so that gamma values are increasing
  reverse(gammaDat.begin(), gammaDat.end());
  reverse(volDat.begin(), volDat.end());
  reverse(energyDat.begin(), energyDat.end());
  reverse(radData.begin(), radData.end());
  return;
}

// Interpolation of area
// Note: This requires xDat to be sorted in descending order.
vec GetInterpolatedValue(double_t x, vec &xDat, vec &y, vec &z, vec &w) {
  int size = xDat.size();
  int i = 0;
  if (x >= xDat[size - 2]) {
    i = size - 2;
  } else {
    while (x > xDat[i + 1])
      i++;
  }
  double_t xL = xDat[i], yL = y[i], zL = z[i], wL = w[i];
  double_t xR = xDat[i + 1], yR = y[i + 1], zR = z[i + 1], wR = w[i + 1];
  if (x < xL) {
    yR = yL;
    zR = zL;
    wR = wL;
  }
  if (x > xR) {
    yL = yR;
    zL = zR;
    wL = wR;
  }
  double_t dydx = (yR - yL) / (xR - xL);
  double_t dzdx = (zR - zL) / (xR - xL);
  double_t dwdx = (wR - wL) / (xR - xL);
  std::vector<double_t> output;
  output.push_back(yL + dydx * (x - xL));
  output.push_back(zL + dzdx * (x - xL));
  output.push_back(wL + dwdx * (x - xL));
  return output;
}

} // namespace OPS
