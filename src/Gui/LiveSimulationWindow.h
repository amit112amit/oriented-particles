#ifndef _LIVESIMULATIONWINDOW_H_
#define _LIVESIMULATIONWINDOW_H_

#include <mutex>
#include <QMainWindow>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QThread>
#include <QVTKOpenGLWidget.h>
#include <vtkOpenGLRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkOpenGLActor.h>
#include <vtkOpenGLProperty.h>
#include <vtkSphereSource.h>
#include <qwt_slider.h>
#include <qwt_scale_engine.h>
#include <qwt_plot.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_curve.h>
#include <qwt_legend.h>
#include <qwt_point_data.h>
#include <qwt_plot_canvas.h>
#include <qwt_text.h>
#include <qwt_symbol.h>
#include <qwt_math.h>
#include <boost/circular_buffer.hpp>
#include "LiveSimulation.h"
#include "ui_LiveSimulationWindow.h"

namespace OPS{

class LiveSimulationWindow : public QMainWindow, public Ui::GUI{
Q_OBJECT
public:
    typedef boost::circular_buffer<double_t> circ_buff;
    explicit LiveSimulationWindow(QWidget *parent=0);
    ~LiveSimulationWindow(){}
signals:
    void startRunning();
    void sceneRefreshed();
    void stopRunning();
    void betaChanged(double);
    void gammaChanged(double);
public slots:
    void refreshVTKSceneAndPlots(int);
    void setUpVTKPipeAndPlotData();
private slots:
    void on__initBtn_clicked();
    void on__startStopBtn_clicked();
    void on__tempSlider_valueChanged(double value);
    void on__gammaSlider_valueChanged(double value);

private:
    vtkSmartPointer<vtkPolyData> _poly;
    LiveSimulation* _worker;
    QThread* _thread;
    circBuffers *_plotDataV;
    std::mutex _mutex;
    QwtPlotGrid *_topGrid, *_bottomGrid;
    QwtPlotCurve *_topCurve, *_bottomCurve;
};

}

#endif // _LIVESIMULATIONWINDOW_H_
