#ifndef _LIVESIMULATION_H_
#define _LIVESIMULATION_H_

#include "ALConstraint.h"
#include "BrownianBody.h"
#include "HelperFunctions.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "Pressure.h"
#include "ViscosityBody.h"
#include <Eigen/Eigenvalues>
#include <QEventLoop>
#include <QObject>
#include <QPointF>
#include <QRectF>
#include <QVector>
#include <algorithm>
#include <boost/circular_buffer.hpp>
#include <limits>
#include <mutex>
#include <qwt_series_data.h>
#include <random>
#include <stdio.h>
#include <string>
#include <vector>
#include <vtkCellData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>

namespace OPS {

// A QwtSeriesData subclass with circular buffer
class QwtCircBuffSeriesData : public QwtSeriesData<QPointF> {
public:
  typedef boost::circular_buffer<qreal> CircBuff;
  // Default constructor
  QwtCircBuffSeriesData(size_t N) {
    _x = new CircBuff(N);
    _y = new CircBuff(N);
    _capacity = N;
    _ymax = _limits.min();
    _ymin = _limits.max();
  }
  ~QwtCircBuffSeriesData() { delete _x, _y; }
  size_t size() const { return _x->size(); }
  QPointF sample(size_t i) const { return QPointF((*_x)[i], (*_y)[i]); }
  QRectF boundingRect() const {
    return QRectF(QPointF(0, _ymax), QPointF(_xmax, _ymin));
  }
  void push_back(qreal x, qreal y) {
    _x->push_back(x);
    _y->push_back(y);
    if (x > _xmax)
      _xmax = x;
    if (y > _ymax)
      _ymax = y;
    if (y < _ymin)
      _ymin = y;
  }
  void clear() {
    _x->clear();
    _y->clear();
    _xmax = 0;
    _ymax = _limits.min();
    _ymin = _limits.max();
  }
  size_t capacity() { return _capacity; }

private:
  CircBuff *_x;
  CircBuff *_y;
  size_t _capacity;
  qreal _ymax = 0, _xmax = 0;
  qreal _ymin = 1e15;
  std::numeric_limits<qreal> _limits;
};

class LiveSimulation : public QObject {
  Q_OBJECT
public:
  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::Vector3d Vector3d;
  typedef std::mt19937 Engine;
  typedef std::normal_distribution<double_t> NormD;
  explicit LiveSimulation() {}
  ~LiveSimulation() { delete _rmsAD, _opsEn, _vol; }
  double_t GetInterpolatedValue(double_t, std::vector<double_t> &);
  void ReadFvKAreaData(std::string);
  double_t GetZeroOpsEn() { return _zeroOpsEnVal; }
  double_t GetZeroRmsAd() { return _zeroRmsAdVal; }
  double_t GetZeroVolume() { return _zeroVolVal; }
  void setInputVTKFile(std::string s) { _inputVTK = s; }
  void setInputZeroData(std::string s) { _inputZeroData = s; }

public slots:
  void Initialize();
  void StartRunning() { _keepRunning = true; }
  void StopRunning() { _keepRunning = false; }
  void UpdateBeta(double b);
  void UpdateGamma(double g);
  void UpdatePressure(double p);
  vtkSmartPointer<vtkPolyData> GetPolyData();
  QwtCircBuffSeriesData *GetRmsAd() { return _rmsAD; }
  QwtCircBuffSeriesData *GetOpsEn() { return _opsEn; }
  QwtCircBuffSeriesData *GetVolume() { return _vol; }
  void SolveOneStep();
  void Reset();
  void LoadState(QString s);
  void SaveState(QString s);
  void SaveScene(QString s);
  void ComputeVoronoi(bool b) { _computeVoronoi = b; };

signals:
  void simulationReady();
  void stepCompleted(int);
  void finished();
  void resetCompeleted();
  void updateZeroOpsEn(double);
  void updateZeroRmsAd(double);
  void updateZeroVolume(double);
  void updatePlotXAxis(int);

private:
  bool _keepRunning = true;
  bool _computeVoronoi = false;
  size_t _step = 0, _N, _resetStepVal = 0;
  double_t _alpha = 2.5e5;
  double_t _beta = 10;
  double_t _gamma = 0.3;
  double_t _pressure = 0.0;
  double_t _zeroRmsAdVal, _zeroOpsEnVal, _zeroVolVal;
  double_t _f = 0.0;
  double_t _R0 = 0.0;
  std::vector<double_t> _gammaDat, _areaDat, _rmsAdDat, _opsEnDat, _volDat;
  std::string _dataOutputFile;
  std::string _inputVTK = "T7.vtk";
  std::string _inputZeroData = "T7_OPS_Asphericity.dat";
  OPSMesh *_ops;
  ViscosityBody *_visco;
  BrownianBody *_brown;
  PressureBody *_pressureBody;
  ExactAreaConstraint *_constraint;
  Model *_model;
  LBFGSBWrapper *_solver;
  std::ofstream _detailedOP;
  VectorXd _x, _g, _prevX, _initialX;
  std::mutex _mut;
  QwtCircBuffSeriesData *_rmsAD, *_opsEn, *_vol;
};

} // namespace OPS
#endif // _LIVESIMULATION_H_
