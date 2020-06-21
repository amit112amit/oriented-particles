#include "LiveSimulationWindow.h"
#include <QApplication>
#include <QSurface>

int main(int argc, char *argv[]) {
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLStereoWidget::defaultFormat());
  QApplication app(argc, argv);
  OPS::LiveSimulationWindow window;
  window.show();
  app.exec();
}
