#include "LiveSimulationWindow.h"

namespace OPS {

LiveSimulationWindow::LiveSimulationWindow(QWidget *parent) {
  setupUi(this);
  // Settings for FvK slider
  _gammaSlider->setScaleEngine(new QwtLogScaleEngine);
  _gammaSlider->setStepAlignment(false);
  _gammaSlider->setBorderWidth(1);
  _gammaSlider->setScale(1.0, 1.0e6);
  _gammaSlider->setTotalSteps(200);
  _gammaSlider->setPageSteps(10);
  _gammaSlider->setScaleMaxMinor(9);
  // Settings for the Pressure Slider
  //_pressureSlider->setScaleEngine( new QwtLogScaleEngine );
  //_pressureSlider->setStepAlignment( false );
  //_pressureSlider->setBorderWidth( 1 );
  //_pressureSlider->setScale( 1.0, 1.0e5 );
  //_pressureSlider->setTotalSteps( 200 );
  //_pressureSlider->setPageSteps( 10 );
  //_pressureSlider->setScaleMaxMinor( 9 );
  // Settings for temperature slider
  _tempSlider->setTotalSteps(200);
  _tempSlider->setValue(0.10);
  _tempValLbl->setText(QString::number(0.10));
  // Create VTK actors and mappers
  _actor = vtkSmartPointer<vtkOpenGLActor>::New();
  _mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  _gActor = vtkSmartPointer<vtkOpenGLActor>::New();
  _gMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  // Create a lookup table to map cell data to colors
  _lookup = vtkSmartPointer<vtkLookupTable>::New();
  _lookup->SetNumberOfTableValues(10);
  _lookup->Build();

  // Fill in a few known colors, the rest will be generated if needed
  _lookup->SetTableValue(0, 0, 0, 0, 1);                // Black
  _lookup->SetTableValue(1, 0.8900, 0.8100, 0.3400, 1); // Banana
  _lookup->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
  _lookup->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
  _lookup->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
  _lookup->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
  _lookup->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
  _lookup->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
  _lookup->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
  _lookup->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock

  // Set up the angle deficit plot
  //_angleDeficitPlot->setTitle("RMS Angle Deficit");
  _angleDeficitPlot->setCanvasBackground(Qt::white);
  _xmax = 5000;
  _ymin1 = 0.10, _ymax1 = 0.5;
  _angleDeficitPlot->setAxisScale(QwtPlot::xBottom, 0, _xmax);
  _angleDeficitPlot->setAxisScale(QwtPlot::yLeft, _ymin1, _ymax1);
  // Create a grid
  _adGrid = new QwtPlotGrid();
  _adGrid->setPen(Qt::black, 1.0, Qt::PenStyle::DashLine);
  // Create a curve
  _angleDeficitCurve = new QwtPlotCurve();
  _angleDeficitCurve->setPen(Qt::blue);
  _angleDeficitCurve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
  // Create a marker
  _adMarker = new QwtPlotMarker();
  _adMarker->setLineStyle(QwtPlotMarker::HLine);
  _adMarker->setLinePen(Qt::red, 1.0, Qt::PenStyle::SolidLine);
  // Attach grid, curve and marker to the plot
  _adGrid->attach(_angleDeficitPlot);
  _angleDeficitCurve->attach(_angleDeficitPlot);
  _adMarker->attach(_angleDeficitPlot);

  // Set up the volume plot
  //_volumePlot->setTitle("Volume");
  _volumePlot->setCanvasBackground(Qt::white);
  _xmax = 5000;
  _ymin3 = 30, _ymax3 = 50;
  _volumePlot->setAxisScale(QwtPlot::xBottom, 0, _xmax);
  _volumePlot->setAxisScale(QwtPlot::yLeft, _ymin3, _ymax3);
  // Create a grid
  _volGrid = new QwtPlotGrid();
  _volGrid->setPen(Qt::black, 1.0, Qt::PenStyle::DashLine);
  // Create a curve
  _volCurve = new QwtPlotCurve();
  _volCurve->setPen(Qt::blue);
  _volCurve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
  // Create a marker
  _volMarker = new QwtPlotMarker();
  _volMarker->setLineStyle(QwtPlotMarker::HLine);
  _volMarker->setLinePen(Qt::red, 1.0, Qt::PenStyle::SolidLine);
  // Attach grid, curve and marker to the plot
  _volGrid->attach(_volumePlot);
  _volCurve->attach(_volumePlot);
  _volMarker->attach(_volumePlot);

  // Set up the energy plot
  //_energyPlot->setTitle("OPS Energy");
  _energyPlot->setCanvasBackground(Qt::white);
  _ymin2 = -500, _ymax2 = 4500;
  _energyPlot->setAxisScale(QwtPlot::xBottom, 0, _xmax);
  _energyPlot->setAxisScale(QwtPlot::yLeft, _ymin2, _ymax2);
  // Create a grid
  _enGrid = new QwtPlotGrid();
  _enGrid->setPen(Qt::black, 1.0, Qt::PenStyle::DashLine);
  // Create a Plot
  _energyCurve = new QwtPlotCurve();
  _energyCurve->setPen(Qt::blue);
  _energyCurve->setRenderHint(QwtPlotItem::RenderAntialiased, true);
  // Create a marker
  _enMarker = new QwtPlotMarker();
  _enMarker->setLineStyle(QwtPlotMarker::HLine);
  _enMarker->setLinePen(Qt::red, 1.0, Qt::PenStyle::SolidLine);
  // Attache grid, curve and marker to the plot
  _enGrid->attach(_energyPlot);
  _energyCurve->attach(_energyPlot);
  _enMarker->attach(_energyPlot);

  // Set up QVTKOpenGLWidget
  _poly = vtkSmartPointer<vtkPolyData>::New();
  _glyphSource = vtkSmartPointer<vtkSphereSource>::New();
  _glyphSource->SetRadius(0.15);
  _glyphSource->SetThetaResolution(10);
  _glyphSource->SetPhiResolution(10);
  auto renderWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  _qVTK->SetRenderWindow(renderWin);

  // Set up connections for worker _thread and object
  _thread = new QThread;
  _worker = new LiveSimulation();
  _worker->moveToThread(_thread);
  connect(_thread, SIGNAL(started()), _worker, SLOT(Initialize()));
  connect(_worker, SIGNAL(finished()), _thread, SLOT(quit()));
  connect(_worker, SIGNAL(finished()), _worker, SLOT(deleteLater()));
  connect(_worker, SIGNAL(simulationReady()), this,
          SLOT(setUpVTKPipeAndPlotData()));
  connect(_worker, SIGNAL(stepCompleted(int)), this,
          SLOT(refreshVTKSceneAndPlots(int)));
  connect(_worker, SIGNAL(resetCompeleted()), this,
          SLOT(resetVTKSceneAndPlots()));
  connect(_worker, SIGNAL(updateZeroOpsEn(double)), this,
          SLOT(updateZeroOpsEnMarker(double)));
  connect(_worker, SIGNAL(updateZeroRmsAd(double)), this,
          SLOT(updateZeroRmsAdMarker(double)));
  connect(_worker, SIGNAL(updateZeroVolume(double)), this,
          SLOT(updateZeroVolumeMarker(double)));
  connect(_worker, SIGNAL(updatePlotXAxis(int)), this,
          SLOT(updatePlotXAxis(int)));
  connect(this, SIGNAL(sceneRefreshed()), _worker, SLOT(SolveOneStep()));
  connect(this, SIGNAL(resetRequested()), _worker, SLOT(Reset()));
  connect(_thread, SIGNAL(finished()), _thread, SLOT(deleteLater()));
  connect(this, SIGNAL(startRunning()), _worker, SLOT(StartRunning()));
  connect(this, SIGNAL(startRunning()), _worker, SLOT(SolveOneStep()));
  connect(this, SIGNAL(stopRunning()), _worker, SLOT(StopRunning()));
  connect(this, SIGNAL(betaChanged(double)), _worker, SLOT(UpdateBeta(double)));
  connect(this, SIGNAL(gammaChanged(double)), _worker,
          SLOT(UpdateGamma(double)));
  connect(this, SIGNAL(pressureChanged(double)), _worker,
          SLOT(UpdatePressure(double)));
  connect(this, SIGNAL(loadStateFile(QString)), _worker,
          SLOT(LoadState(QString)));
  connect(this, SIGNAL(saveStateFile(QString)), _worker,
          SLOT(SaveState(QString)));
  connect(this, SIGNAL(saveVTKFile(QString)), _worker,
          SLOT(SaveScene(QString)));
  connect(this, SIGNAL(showVoronoi(bool)), _worker, SLOT(ComputeVoronoi(bool)));
}

void LiveSimulationWindow::on__initBtn_clicked() {
  if (!_isInitialized) {
    QString inputVTK = QFileDialog::getOpenFileName(
        this, tr("Input VTK File"), "", tr("Legacy VTK (*.vtk)"));
    if (inputVTK.isEmpty())
      return;
    QString inputZeroData = QFileDialog::getOpenFileName(
        this, tr("Input Zero Temperature Data File"), "",
        tr("DAT File(*.dat)"));
    if (inputZeroData.isEmpty())
      return;
    _worker->setInputVTKFile(inputVTK.toStdString());
    _worker->setInputZeroData(inputZeroData.toStdString());
    _initBtn->setText(QString("Reset"));
    _isInitialized = true;
    if (_thread)
      _thread->start();
  } else {
    emit resetRequested();
  }
}

void LiveSimulationWindow::updatePlotXAxis(int i) {
  _xmin = i;
  _xmax = _xmin + 5000;
  _angleDeficitPlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
  _energyPlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
  _volumePlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
  _angleDeficitPlot->replot();
  _energyPlot->replot();
  _volumePlot->replot();
}

void LiveSimulationWindow::setUpVTKPipeAndPlotData() {
  // Set up QVTKOpenGLRenderer
  // std::lock_guard<std::mutex> lock(_mutex);
  _poly->DeepCopy(_worker->GetPolyData());
  _poly->GetCellData()->SetActiveScalars(_colorarray.c_str());
  auto glyph = vtkSmartPointer<vtkGlyph3D>::New();
  auto renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
  // Set the shell mapper and actor
  _mapper->SetInputData(_poly);
  _mapper->SetScalarRange(0, 9);
  _mapper->SetLookupTable(_lookup);
  //_actor->GetProperty()->EdgeVisibilityOn();
  //_actor->GetProperty()->SetEdgeColor(0.3,1,1);
  // Set the glyph mapper and actor
  glyph->SetInputData(_poly);
  glyph->SetSourceConnection(_glyphSource->GetOutputPort());
  _gMapper->SetInputConnection(glyph->GetOutputPort());
  // We need to turn off some actors
  if (_showVoronoi) {
    _gActor->VisibilityOff();
    _mapper->ScalarVisibilityOn();
    _mapper->SetScalarRange(0, 9);
    _mapper->SetLookupTable(_lookup);
  } else {
    _gActor->VisibilityOn();
    _gMapper->ScalarVisibilityOff();
    _gActor->GetProperty()->SetColor(0.6666667, 0., 0.);
    _mapper->ScalarVisibilityOff();
    _actor->GetProperty()->SetColor(0.6666667, 1., 1.);
  }
  // Set mappers
  _actor->SetMapper(_mapper);
  _gActor->SetMapper(_gMapper);
  // Add the two actors to renderer
  renderer->AddViewProp(_actor);
  renderer->AddViewProp(_gActor);
  renderer->SetBackground(0.318, 0.341, 0.431);
  _qVTK->GetRenderWindow()->AddRenderer(renderer);
  _qVTK->GetRenderWindow()->Render();
  // Also get addresses of plotting data circular buffers;
  _rmsAd = _worker->GetRmsAd();
  _totEn = _worker->GetOpsEn();
  _vol = _worker->GetVolume();
  _angleDeficitCurve->setSamples(_rmsAd);
  _volCurve->setSamples(_vol);
  _energyCurve->setSamples(_totEn);
  _adMarker->setYValue(_worker->GetZeroRmsAd());
  _volMarker->setYValue(_worker->GetZeroVolume());
  _enMarker->setYValue(_worker->GetZeroOpsEn());
  _volumePlot->replot();
  _angleDeficitPlot->replot();
  _energyPlot->replot();
}

void LiveSimulationWindow::refreshVTKSceneAndPlots(int s) {
  //{
  // Update the VTK scene
  // std::lock_guard<std::mutex> lock(_mutex);
  _poly->DeepCopy(_worker->GetPolyData());
  _poly->GetCellData()->SetActiveScalars(_colorarray.c_str());
  _poly->Modified();
  if (_showVoronoi) {
    _gActor->VisibilityOff();
    _mapper->ScalarVisibilityOn();
    _mapper->SetScalarRange(0, 9);
    _mapper->SetLookupTable(_lookup);
  } else {
    _gActor->VisibilityOn();
    _gMapper->ScalarVisibilityOff();
    _gActor->GetProperty()->SetColor(0.6666667, 0., 0.);
    _mapper->ScalarVisibilityOff();
    _actor->GetProperty()->SetColor(0.6666667, 1., 1.);
  }
  _qVTK->GetRenderWindow()->Render();
  // Update the plots
  qreal x1, y1, x2, y2;
  _angleDeficitCurve->boundingRect().getCoords(&x1, &y1, &x2, &y2);
  if (x2 > _xmax) {
    _xmin += 1000;
    _xmax += 1000;
    _angleDeficitPlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
    _energyPlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
    _volumePlot->setAxisScale(QwtPlot::xBottom, _xmin, _xmax);
  }
  if (y1 > _ymax1 || y2 < _ymin1) {
    _ymin1 = std::min(0.75 * y2, _ymin1);
    _ymax1 = std::max(1.5 * y1, _ymax1);
    _angleDeficitPlot->setAxisScale(QwtPlot::yLeft, _ymin1, _ymax1);
  }
  _volCurve->boundingRect().getCoords(&x1, &y1, &x2, &y2);
  if (y1 > _ymax3 || y2 < _ymin3) {
    _ymin3 = std::min(0.75 * y2, _ymin3);
    _ymax3 = std::max(1.5 * y1, _ymax3);
    _volumePlot->setAxisScale(QwtPlot::yLeft, _ymin3, _ymax3);
  }
  _energyCurve->boundingRect().getCoords(&x1, &y1, &x2, &y2);
  if (y1 > _ymax2 || y2 < _ymin2) {
    _ymin2 = std::min(0.75 * y2, _ymin2);
    _ymax2 = std::max(1.5 * y1, _ymax2);
    _energyPlot->setAxisScale(QwtPlot::yLeft, _ymin2, _ymax2);
  }
  // Finally update the plot
  _angleDeficitPlot->replot();
  _energyPlot->replot();
  _volumePlot->replot();
  //}
  // Update time step
  _timeStepValLbl->setText(QString::number(s));
  emit sceneRefreshed();
}

void LiveSimulationWindow::resetVTKSceneAndPlots() {
  _poly->DeepCopy(_worker->GetPolyData());
  _poly->GetCellData()->SetActiveScalars(_colorarray.c_str());
  _poly->Modified();
  if (_showVoronoi) {
    _gActor->VisibilityOff();
    _mapper->ScalarVisibilityOn();
    _mapper->SetScalarRange(0, 9);
    _mapper->SetLookupTable(_lookup);
  } else {
    _gActor->VisibilityOn();
    _gMapper->ScalarVisibilityOff();
    _gActor->GetProperty()->SetColor(0.6666667, 0., 0.);
    _mapper->ScalarVisibilityOff();
    _actor->GetProperty()->SetColor(0.6666667, 1., 1.);
  }
  _qVTK->GetRenderWindow()->Render();
  _timeStepValLbl->setText(QString::number(int(_xmin)));
  _startStopBtn->setText(QString("Start"));
}

void LiveSimulationWindow::on__startStopBtn_clicked() {
  QString txt = "Start";
  if (QString::compare(_startStopBtn->text(), txt) == 0) {
    _startStopBtn->setText(QString("Stop"));
    _rmsAd = _worker->GetRmsAd();
    _totEn = _worker->GetOpsEn();
    _vol = _worker->GetVolume();
    emit startRunning();
  } else {
    _startStopBtn->setText(QString("Start"));
    emit stopRunning();
  }
}

void LiveSimulationWindow::on__tempSlider_valueChanged(double value) {
  emit betaChanged(value);
}

void LiveSimulationWindow::on__gammaSlider_valueChanged(double value) {
  emit gammaChanged(value);
}

void LiveSimulationWindow::on__pressureSlider_valueChanged(double value) {
  emit pressureChanged(value);
}

void LiveSimulationWindow::updateZeroOpsEnMarker(double o) {
  _enMarker->setYValue(o);
  _energyPlot->replot();
}

void LiveSimulationWindow::updateZeroRmsAdMarker(double r) {
  _adMarker->setYValue(r);
  _angleDeficitPlot->replot();
}

void LiveSimulationWindow::updateZeroVolumeMarker(double r) {
  _volMarker->setYValue(r);
  _volumePlot->replot();
}

void LiveSimulationWindow::on_pressureCheckBox_clicked(bool checked) {
  if (checked) {
    _pressureSlider->setEnabled(true);
    _pressureValLbl->setText(QString::number(_pressureSlider->value()));
    emit pressureChanged(_pressureSlider->value());
  } else {
    _pressureSlider->setEnabled(false);
    _pressureValLbl->setText(QString::number(0.0));
    emit pressureChanged(0.0);
  }
}

void LiveSimulationWindow::on_voronoiCheckBox_clicked(bool checked) {
  _showVoronoi = checked;
  emit showVoronoi(checked);
}

void LiveSimulationWindow::on_actionLoadState_triggered() {
  if (!_isInitialized)
    on__initBtn_clicked();
  if (_isInitialized) {
    QString fileName = QFileDialog::getOpenFileName(
        this, tr("Open Simulation State File"), "", tr("DAT File (*.dat)"));
    if (!fileName.isEmpty()) {
      SimulationState state =
          SimulationState::readFromFile(fileName.toStdString());
      double_t g = state.getGamma();
      double_t t = 1.0 / state.getBeta();
      _gammaSlider->setValue(g);
      _fvkValLbl->setText(QString::number(g));
      _tempSlider->setValue(t);
      _tempValLbl->setText(QString::number(t));
      _timeStepValLbl->setText(QString::number(state.getStep() + 1));
      emit loadStateFile(fileName);
    }
  }
}

void LiveSimulationWindow::on_actionSaveState_triggered() {
  if (_isInitialized) {
    _startStopBtn->setText(QString("Start"));
    emit stopRunning();
    QString saveFileName = QFileDialog::getSaveFileName(
        this, tr("Save Simulation State"), "", tr("DAT File (*.dat)"));
    if (!saveFileName.isEmpty()) {
      QString ext = ".dat";
      if (!saveFileName.endsWith(ext))
        saveFileName.append(ext);
      emit saveStateFile(saveFileName);
    }
  } else {
    QMessageBox messageBox;
    messageBox.critical(0, "Error", "Simulation has not been started!");
    messageBox.setFixedSize(500, 200);
  }
}

void LiveSimulationWindow::on_actionExportScene_triggered() {
  QString saveFileName = QFileDialog::getSaveFileName(
      this, tr("Export to VTK file"), "", tr("VTK File (*.vtk)"));
  if (!saveFileName.isEmpty()) {
    QString ext = ".vtk";
    if (!saveFileName.endsWith(ext))
      saveFileName.append(ext);
    emit saveVTKFile(saveFileName);
  }
}

} // namespace OPS
