#ifndef _LIVESIMULATION_H_
#define _LIVESIMULATION_H_

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
#include "ALConstraint.h"
#include "BrownianBody.h"
#include "LBFGSBWrapper.h"
#include "Model.h"
#include "OPSMesh.h"
#include "ViscosityBody.h"

namespace OPS{

struct circBuffers{
    boost::circular_buffer<double_t>* x;
    boost::circular_buffer<double_t>* a;
    boost::circular_buffer<double_t>* b;
    double_t* getA(){return &(*a)[0];}
    double_t* getB(){return &(*b)[0];}
    double_t* getX(){return &(*x)[0];}
};

class LiveSimulation: public QObject{
    Q_OBJECT
public:
    typedef boost::circular_buffer<double_t> circ_buff;
    typedef Eigen::VectorXd VectorXd;
    typedef Eigen::Vector3d Vector3d;
    explicit LiveSimulation(){}
    ~LiveSimulation(){}
    double_t getInterpolatedArea(double_t);
    void readFvKAreaData(std::string);
public slots:
    void Initialize();
    void StartRunning(){_keepRunning = true;}
    void StopRunning(){_keepRunning = false;}
    void updateBeta(double b);
    void updateGamma(double g);
    vtkSmartPointer<vtkPolyData> getPolyData();
    circBuffers* getPlotData(){return _plotData;}
    void SolveOneStep();

signals:
    void simulationReady();
    void stepCompleted(int);
    void finished();
private:
    bool _keepRunning = true;
    size_t _step = 0, _N;
    double_t _alpha = 2.5e5;
    double_t _beta = 10;
    double_t _gamma = 0.01;
    double_t _f = 0.0;
    std::vector<double_t> _gammaDat, _areaDat;
    OPSMesh* _ops;
    ViscosityBody* _visco;
    BrownianBody* _brown;
    ExactAreaConstraint* _constraint;
    Model* _model;
    LBFGSBWrapper* _solver;
    ofstream _detailedOP;
    VectorXd _x, _g, _prevX, _initialX;
    std::mutex _mut;
    circBuffers *_plotData;
};

}
#endif // _LIVESIMULATION_H_
