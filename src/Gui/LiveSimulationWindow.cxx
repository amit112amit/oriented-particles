#include "LiveSimulationWindow.h"

namespace OPS{

LiveSimulationWindow::LiveSimulationWindow(QWidget *parent){
    setupUi(this);
    // Settings for FvK slider
    _gammaSlider->setOrientation(Qt::Horizontal);
    _gammaSlider->setTrough( true );
    _gammaSlider->setScaleEngine( new QwtLogScaleEngine );
    _gammaSlider->setStepAlignment( false );
    _gammaSlider->setBorderWidth( 1 );
    _gammaSlider->setScale( 1.0e-2, 1.0e4 );
    _gammaSlider->setTotalSteps( 200 );
    _gammaSlider->setPageSteps( 10 );
    _gammaSlider->setScaleMaxMinor( 9 );
    // Settings for temperature slider
    _tempSlider->setTotalSteps( 200 );
    _tempSlider->setValue(0.01);
    _tempValLbl->setText(QString::number(0.01));

    // Set up the plots
    _topPlot->setTitle("RMS Angle Deficit");
    _topPlot->setCanvasBackground(Qt::white);
    _topPlot->setAxisScale(QwtPlot::yLeft,0.17,0.3);
    _topPlot->axisAutoScale(QwtPlot::yLeft);
    _topPlot->setAxisScale(QwtPlot::xBottom,0,10000);
    _topCurve = new QwtPlotCurve();
    _topGrid = new QwtPlotGrid();
    _topCurve->setPen(Qt::blue);
    _topCurve->setRenderHint(QwtPlotItem::RenderAntialiased,true);
    _topCurve->attach(_topPlot);
    _topGrid->attach(_topPlot);

    _bottomPlot->setTitle("OPS Energy");
    _bottomPlot->setCanvasBackground(Qt::white);
    _bottomPlot->setAxisScale(QwtPlot::yLeft,200,500);
    _bottomPlot->axisAutoScale(QwtPlot::yLeft);
    _bottomPlot->setAxisScale(QwtPlot::xBottom,0,10000);
    _bottomCurve = new QwtPlotCurve();
    _bottomGrid = new QwtPlotGrid();
    _bottomCurve->setPen(Qt::blue);
    _bottomCurve->setRenderHint(QwtPlotItem::RenderAntialiased,true);
    _bottomCurve->attach(_bottomPlot);
    _bottomGrid->attach(_bottomPlot);

    // Set up QVTKOpenGLWidget
    auto renderWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    _qVTK->SetRenderWindow(renderWin);

    //Set up connections for worker _thread and object
    _thread = new QThread;
    _worker = new LiveSimulation();
    _worker->moveToThread(_thread);
    connect(_thread, SIGNAL(started()), _worker, SLOT(Initialize()));
    connect(_worker, SIGNAL(finished()), _thread, SLOT(quit()));
    connect(_worker, SIGNAL(finished()), _worker, SLOT(deleteLater()));
    connect(_worker, SIGNAL(simulationReady()), this, SLOT(setUpVTKPipeAndPlotData()));
    connect(_worker, SIGNAL(stepCompleted(int)), this, SLOT(refreshVTKSceneAndPlots(int)));
    connect(this, SIGNAL(sceneRefreshed()), _worker, SLOT(SolveOneStep()));
    connect(_thread, SIGNAL(finished()), _thread, SLOT(deleteLater()));
    connect(this, SIGNAL(startRunning()), _worker, SLOT(StartRunning()));
    connect(this, SIGNAL(startRunning()), _worker, SLOT(SolveOneStep()));
    connect(this, SIGNAL(stopRunning()), _worker, SLOT(StopRunning()));
    connect(this, SIGNAL(betaChanged(double)), _worker, SLOT(updateBeta(double)));
    connect(this, SIGNAL(gammaChanged(double)), _worker, SLOT(updateGamma(double)));
}

void LiveSimulationWindow::on__initBtn_clicked(){
    if(_thread)
        _thread->start();
}

void LiveSimulationWindow::setUpVTKPipeAndPlotData(){
    // Set up QVTKOpenGLRenderer
    //std::lock_guard<std::mutex> lock(_mutex);
    _poly = _worker->getPolyData();
    auto _mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    _mapper->SetInputData(_poly);
    _mapper->ScalarVisibilityOff();
    auto actor = vtkSmartPointer<vtkOpenGLActor>::New();
    actor->SetMapper(_mapper);
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetEdgeColor(0,0,0.5);
    auto renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
    renderer->AddViewProp(actor);
    renderer->SetBackground(0.318,0.341,0.431);
    _qVTK->GetRenderWindow()->AddRenderer(renderer);
    _qVTK->GetRenderWindow()->Render();
    // Also get addresses of plotting data circular buffers;
    _plotDataV = _worker->getPlotData();
}

void LiveSimulationWindow::refreshVTKSceneAndPlots(int s){
    //{
        // Update the VTK scene
        //std::lock_guard<std::mutex> lock(_mutex);
        _poly->Modified();
        _qVTK->GetRenderWindow()->Render();
        // Update the plots
        if(s > 10000){
            _topPlot->axisAutoScale(QwtPlot::xBottom);
            _bottomPlot->axisAutoScale(QwtPlot::xBottom);
        }
        _topCurve->setRawSamples(_plotDataV->getX(),
                                 _plotDataV->getA(),s);
        _topPlot->replot();
        _bottomCurve->setRawSamples(_plotDataV->getX(),
                                    _plotDataV->getB(),s);
        _bottomPlot->replot();
    //}
    _timeStepValLbl->setText(QString::number(s));
    emit sceneRefreshed();
}

void LiveSimulationWindow::on__startStopBtn_clicked(){
    QString txt = "Start";
    if(QString::compare(_startStopBtn->text(),txt) == 0){
        _startStopBtn->setText(QString("Stop"));
        emit startRunning();
    }
    else{
        _startStopBtn->setText(QString("Start"));
        emit stopRunning();
    }
}

void LiveSimulationWindow::on__tempSlider_valueChanged(double value){
    emit betaChanged(value);
}

void LiveSimulationWindow::on__gammaSlider_valueChanged(double value){
    emit gammaChanged(value);
}

}
