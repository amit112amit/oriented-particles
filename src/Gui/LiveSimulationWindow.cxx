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
    _tempSlider->setValue(0.10);
    _tempValLbl->setText(QString::number(0.10));

    // Set up the plots

    _topPlot->setTitle("RMS Angle Deficit");
    _topPlot->setCanvasBackground(Qt::white);
    _xmax = 5000; _ymin1 = 0.10, _ymax1 = 0.5;
    _topPlot->setAxisScale(QwtPlot::xBottom,0,_xmax);
    _topPlot->setAxisScale(QwtPlot::yLeft,_ymin1,_ymax1);
    // Create a grid
    _topGrid = new QwtPlotGrid();
    _topGrid->setPen(Qt::black,1.0,Qt::PenStyle::DashLine);
    // Create a curve
    _topCurve = new QwtPlotCurve();
    _topCurve->setPen(Qt::blue);
    _topCurve->setRenderHint(QwtPlotItem::RenderAntialiased,true);
    // Create a marker
    _topMarker = new QwtPlotMarker();
    _topMarker->setLineStyle(QwtPlotMarker::HLine);
    _topMarker->setLinePen(Qt::red,1.0,Qt::PenStyle::SolidLine);
    // Attach grid, curve and marker to the plot
    _topGrid->attach(_topPlot);
    _topCurve->attach(_topPlot);
    _topMarker->attach(_topPlot);

    _bottomPlot->setTitle("OPS Energy");
    _bottomPlot->setCanvasBackground(Qt::white);
    _ymin2 = -500, _ymax2 = 4500;
    _bottomPlot->setAxisScale(QwtPlot::xBottom,0,_xmax);
    _bottomPlot->setAxisScale(QwtPlot::yLeft,_ymin2,_ymax2);
    // Create a grid
    _bottomGrid = new QwtPlotGrid();
    _bottomGrid->setPen(Qt::black,1.0,Qt::PenStyle::DashLine);
    // Create a Plot
    _bottomCurve = new QwtPlotCurve();
    _bottomCurve->setPen(Qt::blue);
    _bottomCurve->setRenderHint(QwtPlotItem::RenderAntialiased,true);
    // Create a marker
    _bottomMarker = new QwtPlotMarker();
    _bottomMarker->setLineStyle(QwtPlotMarker::HLine);
    _bottomMarker->setLinePen(Qt::red,1.0,Qt::PenStyle::SolidLine);
    // Attache grid, curve and marker to the plot
    _bottomGrid->attach(_bottomPlot);
    _bottomCurve->attach(_bottomPlot);
    _bottomMarker->attach(_bottomPlot);

    // Set up QVTKOpenGLWidget
    _poly = vtkSmartPointer<vtkPolyData>::New();
    _glyphSource = vtkSmartPointer<vtkSphereSource>::New();
    _glyphSource->SetRadius(0.15);
    _glyphSource->SetThetaResolution(10);
    _glyphSource->SetPhiResolution(10);
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
    connect(_worker, SIGNAL(resetCompeleted()), this, SLOT(resetVTKSceneAndPlots()));
    connect(_worker, SIGNAL(updateZeroOpsEn(double)), this, SLOT(updateZeroOpsEnMarker(double)));
    connect(_worker, SIGNAL(updateZeroRmsAd(double)), this, SLOT(updateZeroRmsAdMarker(double)));
    connect(this, SIGNAL(sceneRefreshed()), _worker, SLOT(SolveOneStep()));
    connect(this, SIGNAL(resetRequested()), _worker, SLOT(Reset()));
    connect(_thread, SIGNAL(finished()), _thread, SLOT(deleteLater()));
    connect(this, SIGNAL(startRunning()), _worker, SLOT(StartRunning()));
    connect(this, SIGNAL(startRunning()), _worker, SLOT(SolveOneStep()));
    connect(this, SIGNAL(stopRunning()), _worker, SLOT(StopRunning()));
    connect(this, SIGNAL(betaChanged(double)), _worker, SLOT(UpdateBeta(double)));
    connect(this, SIGNAL(gammaChanged(double)), _worker, SLOT(UpdateGamma(double)));
}

void LiveSimulationWindow::on__initBtn_clicked(){
    QString txt = "Initialize";
    if(QString::compare(_initBtn->text(),txt) == 0){
        _initBtn->setText(QString("Reset"));
        if(_thread)
            _thread->start();
    }
    else{
        emit resetRequested();
    }
}

void LiveSimulationWindow::setUpVTKPipeAndPlotData(){
    // Set up QVTKOpenGLRenderer
    //std::lock_guard<std::mutex> lock(_mutex);
    _poly->DeepCopy(_worker->GetPolyData());
    auto glyph = vtkSmartPointer<vtkGlyph3D>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    auto gMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    auto actor = vtkSmartPointer<vtkOpenGLActor>::New();
    auto gActor = vtkSmartPointer<vtkOpenGLActor>::New();
    auto renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
    // Set the shell mapper and actor
    mapper->SetInputData(_poly);
    mapper->ScalarVisibilityOff();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0.6666667,1.,1.);
    //actor->GetProperty()->EdgeVisibilityOn();
    //actor->GetProperty()->SetEdgeColor(0.3,1,1);
    // Set the glyph mapper and actor
    glyph->SetInputData(_poly);
    glyph->SetSourceConnection(_glyphSource->GetOutputPort());
    gMapper->SetInputConnection(glyph->GetOutputPort());
    gMapper->ScalarVisibilityOff();
    gActor->SetMapper(gMapper);
    gActor->GetProperty()->SetColor(0.6666667,0.,0.);
    // Add the two actors to renderer
    renderer->AddViewProp(actor);
    renderer->AddViewProp(gActor);
    renderer->SetBackground(0.318,0.341,0.431);
    _qVTK->GetRenderWindow()->AddRenderer(renderer);
    _qVTK->GetRenderWindow()->Render();
    // Also get addresses of plotting data circular buffers;
    _rmsAd = _worker->GetRmsAd();
    _opsEn = _worker->GetOpsEn();
    _topCurve->setSamples(_rmsAd);
    _bottomCurve->setSamples(_opsEn);
    _topMarker->setYValue(_worker->GetZeroRmsAd());
    _bottomMarker->setYValue(_worker->GetZeroOpsEn());
}

void LiveSimulationWindow::refreshVTKSceneAndPlots(int s){
    //{
    // Update the VTK scene
    //std::lock_guard<std::mutex> lock(_mutex);
    _poly->DeepCopy(_worker->GetPolyData());
    _poly->Modified();
    _qVTK->GetRenderWindow()->Render();
    // Update the plots
    qreal x1,y1,x2,y2;
    _topCurve->boundingRect().getCoords(&x1,&y1,&x2,&y2);
    if(x2 > _xmax){
        _xmin += 1000;
        _xmax += 1000;
        _topPlot->setAxisScale(QwtPlot::xBottom,_xmin,_xmax);
        _bottomPlot->setAxisScale(QwtPlot::xBottom,_xmin,_xmax);
    }
    if(y1 > _ymax1 || y2 < _ymin1){
        _ymin1 = std::min(1.5*y2,_ymin1);
        _ymax1 = std::max(1.5*y1,_ymax1);
        _topPlot->setAxisScale(QwtPlot::yLeft,_ymin1,_ymax1);
    }
    _bottomCurve->boundingRect().getCoords(&x1,&y1,&x2,&y2);
    if(y1 > _ymax2 || y2 < _ymin2){
        _ymin2 = std::min(1.5*y2,_ymin2);
        _ymax2 = std::max(1.5*y1,_ymax2);
        _bottomPlot->setAxisScale(QwtPlot::yLeft,_ymin2,_ymax2);
    }
    // Finally update the plot
    _topPlot->replot();
    _bottomPlot->replot();
    //}
    // Update time step
    _timeStepValLbl->setText(QString::number(s));
    emit sceneRefreshed();
}

void LiveSimulationWindow::resetVTKSceneAndPlots(){
    _poly->DeepCopy(_worker->GetPolyData());
    _poly->Modified();
    _qVTK->GetRenderWindow()->Render();
    _timeStepValLbl->setText(QString::number(0));
    _startStopBtn->setText(QString("Start"));
    _xmin = 0;
    _xmax = 5000;
    _topPlot->setAxisScale(QwtPlot::xBottom,_xmin,_xmax);
    _bottomPlot->setAxisScale(QwtPlot::xBottom,_xmin,_xmax);
    _topPlot->replot();
    _bottomPlot->replot();
}

void LiveSimulationWindow::on__startStopBtn_clicked(){
    QString txt = "Start";
    if(QString::compare(_startStopBtn->text(),txt) == 0){
        _startStopBtn->setText(QString("Stop"));
        _rmsAd = _worker->GetRmsAd();
        _opsEn = _worker->GetOpsEn();
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

void LiveSimulationWindow::updateZeroOpsEnMarker(double o){
    _bottomMarker->setYValue(o);
    _bottomPlot->replot();
}

void LiveSimulationWindow::updateZeroRmsAdMarker(double r){
    _topMarker->setYValue(r);
    _topPlot->replot();
}

}
