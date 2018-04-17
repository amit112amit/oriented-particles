#include <QApplication>
#include <QSurface>
#include "LiveSimulationWindow.h"

int main(int argc, char* argv[]){
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
    QApplication app(argc,argv);
    OPS::LiveSimulationWindow window;
    window.show();
    app.exec();
}
