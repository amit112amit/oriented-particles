#ifndef _ROUGHNESSWIDGET_H
#define _ROUGHNESSWIDGET_H

#include "ui_RoughnessGUI.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkOpenGLActor.h"
#include "vtkOpenGLProperty.h"
#include "vtkOpenGLRenderer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkRendererCollection.h"
#include "vtkSmartPointer.h"
#include <Eigen/Dense>
#include <QMainWindow>
//#include <QVTKOpenGLWidget.h>
#include <QVTKOpenGLNativeWidget.h>
#include <cstdlib>
#include <math.h>

class RoughnessWidget : public QMainWindow, public Ui::GUI {
  Q_OBJECT
public:
  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::Vector3d Vector3d;
  explicit RoughnessWidget(QWidget *parent = 0);
  ~RoughnessWidget() {}
signals:
  void updateScene1();
  void updateScene2();
private slots:
  double_t calculateAsphericity(vtkSmartPointer<vtkPolyData> p);
  double_t calculateRMSAd(vtkSmartPointer<vtkPolyData> p);
  void on_slider1_valueChanged(double position);
  void on_slider2_valueChanged(double position);

private:
  vtkSmartPointer<vtkPolyData> poly1;
  vtkSmartPointer<vtkPolyData> poly2;
  vtkSmartPointer<vtkPolyData> poly1_0;
  vtkSmartPointer<vtkPolyData> poly2_0;
};

#endif // _ROUGHNESSWIDGET_H
