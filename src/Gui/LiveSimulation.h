#ifndef _LIVESIMULATION_H_
#define _LIVESIMULATION_H_

#include <limits>
#include <mutex>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <boost/circular_buffer.hpp>
#include <Eigen/Eigenvalues>
#include <QObject>
#include <QEventLoop>
#include <QVector>
#include <QPointF>
#include <QRectF>
#include <qwt_series_data.h>
#include "ALConstraint.h"
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"
#include "Pressure.h"

namespace OPS{

// A QwtSeriesData subclass with circular buffer
class QwtCircBuffSeriesData: public QwtSeriesData<QPointF>{
public:
    QwtCircBuffSeriesData(size_t N){
        _x = new boost::circular_buffer<qreal>(N);
        _y = new boost::circular_buffer<qreal>(N);
        _capacity = N;
        _ymax = _limits.min();
        _ymin = _limits.max();
    }
    ~QwtCircBuffSeriesData(){
        delete _x, _y;
    }
    size_t size() const{
        return _x->size();
    }
    QPointF sample(size_t i) const{
        return QPointF((*_x)[i],(*_y)[i]);
    }
    QRectF boundingRect() const{
        return QRectF(QPointF(0,_ymax),QPointF(_xmax,_ymin));
    }
    void push_back(qreal x, qreal y){
        _x->push_back(x);
        _y->push_back(y);
        if(x > _xmax)
            _xmax = x;
        if(y > _ymax)
            _ymax = y;
        if(y < _ymin)
            _ymin = y;
    }
    void clear(){
        _x->clear();
        _y->clear();
        _xmax = 0;
        _ymax = _limits.min();
        _ymin = _limits.max();
    }
    size_t capacity(){
        return _capacity;
    }
private:
    boost::circular_buffer<qreal> *_x;
    boost::circular_buffer<qreal> *_y;
    size_t _capacity;
    qreal _ymax = 0, _xmax = 0;
    qreal _ymin = 1e15;
    std::numeric_limits<qreal> _limits;
};

class LiveSimulation: public QObject{
    Q_OBJECT
public:
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Vector3d Vector3d;
    explicit LiveSimulation(){}
    ~LiveSimulation(){}
    double_t GetInterpolatedValue(double_t, std::vector<double_t>&);
    void ReadFvKAreaData(std::string);
    double_t GetZeroOpsEn(){return _zeroOpsEnVal;}
    double_t GetZeroRmsAd(){return _zeroRmsAdVal;}
    double_t GetZeroVolume(){return _zeroVolVal;}
    void setInputVTKFile(std::string s){_inputVTK = s;}
    void setInputZeroData(std::string s){_inputZeroData = s;}

public slots:
    void Initialize();
    void StartRunning(){_keepRunning = true;}
    void StopRunning(){_keepRunning = false;}
    void UpdateBeta(double b);
    void UpdateGamma(double g);
    void UpdatePressure(double p);
    vtkSmartPointer<vtkPolyData> GetPolyData();
    QwtCircBuffSeriesData* GetRmsAd(){return _rmsAD;}
    QwtCircBuffSeriesData* GetOpsEn(){return _opsEn;}
    QwtCircBuffSeriesData* GetVolume(){return _vol;}
    void SolveOneStep();
    void Reset();
    void LoadState( QString s );

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
    size_t _step = 0, _N, _resetStepVal = 0;
    double_t _alpha = 2.5e5;
    double_t _beta = 10;
    double_t _gamma = 0.01;
    double_t _pressure = 0.0;
    double_t _zeroRmsAdVal, _zeroOpsEnVal, _zeroVolVal;
    double_t _f = 0.0;
    std::vector<double_t> _gammaDat, _areaDat,
    _rmsAdDat, _opsEnDat, _volDat;
    std::string _dataOutputFile;
    std::string _inputVTK = "T7.vtk";
    std::string _inputZeroData = "T7_OPS_Asphericity.dat";
    OPSMesh* _ops;
    ViscosityBody* _visco;
    BrownianBody* _brown;
    InternalPressure* _pressureBody;
    ExactAreaConstraint* _constraint;
    Model* _model;
    LBFGSBWrapper* _solver;
    ofstream _detailedOP;
    VectorXd _x, _g, _prevX, _initialX;
    std::mutex _mut;
    QwtCircBuffSeriesData *_rmsAD, *_opsEn, *_vol;
};

}
#endif // _LIVESIMULATION_H_
