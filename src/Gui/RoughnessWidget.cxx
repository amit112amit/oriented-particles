#include "RoughnessWidget.h"

RoughnessWidget::RoughnessWidget(QWidget *parent) : QMainWindow(parent){
    setupUi(this);
    srand(time(NULL));

    // put poly1 in one window
    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName("T7-relaxed.vtk");
    reader->Update();
    poly1 = reader->GetOutput();
    poly1_0 = vtkSmartPointer<vtkPolyData>::New();
    poly1_0->DeepCopy(poly1);
    auto mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->SetInputData(poly1);
    mapper1->ScalarVisibilityOff();
    auto actor1 = vtkSmartPointer<vtkOpenGLActor>::New();
    actor1->SetMapper(mapper1);
    actor1->GetProperty()->SetEdgeVisibility(1);
    actor1->GetProperty()->SetColor(1,1,1);
    actor1->GetProperty()->SetEdgeColor(0,0,0.5);
    auto Ren1 = vtkSmartPointer<vtkOpenGLRenderer>::New();
    Ren1->AddActor(actor1);
    Ren1->SetBackground(0.318,0.341,0.431);
    auto renwin1 = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    renwin1->AddRenderer(Ren1);
    qVTK1->SetRenderWindow(renwin1);


    // put poly2 in other window
    auto mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    poly2 = vtkSmartPointer<vtkPolyData>::New();
    poly2->DeepCopy(poly1);
    poly2_0 = vtkSmartPointer<vtkPolyData>::New();
    poly2_0->DeepCopy(poly2);
    mapper2->SetInputData(poly2);
    mapper2->ScalarVisibilityOff();
    auto actor = vtkSmartPointer<vtkOpenGLActor>::New();
    actor->SetMapper(mapper2);
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetColor(1,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0.5);
    auto Ren2 = vtkSmartPointer<vtkOpenGLRenderer>::New();
    Ren2->AddActor(actor);
    Ren2->SetBackground(0.318,0.341,0.431);
    auto renwin2 = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    renwin2->AddRenderer(Ren2);
    qVTK2->SetRenderWindow(renwin2);

    //Set Asphericity and RMS Angle Deficit
    double_t val = calculateAsphericity(poly1);
    QString txt = QString("Asphericity: %1").arg(val,0,'g',3);
    label1->setText(txt);
    val = calculateAsphericity(poly2);
    txt = QString("Asphericity: %1").arg(val,0,'g',3);
    label2->setText(txt);

    val = calculateRMSAd(poly1);
    txt = QString("RMS Angle Deficit: %1").arg(val,0,'g',3);
    label3->setText(txt);
    val = calculateRMSAd(poly2);
    txt = QString("RMS Angle Deficit: %1").arg(val,0,'g',3);
    label4->setText(txt);
}

void RoughnessWidget::on_slider2_valueChanged(double position){
    auto N = poly2->GetNumberOfPoints();
    //for(auto i=N-12; i < N; ++i){
    for(auto i=0; i < N; ++i){
        Vector3d p;
        int randomBit = rand() & 1;
        randomBit = (randomBit > 0)? 1 : -1;
        poly2_0->GetPoint(i,&p(0));
        p -= position*randomBit*p;
        poly2->GetPoints()->SetPoint(i,&p(0));
    }
    double_t val = calculateAsphericity(poly2);
    QString txt = QString("Asphericity: %1").arg(val,0,'g',3);
    label2->setText(txt);
    val = calculateRMSAd(poly2);
    txt = QString("RMS Angle Deficit: %1").arg(val,0,'g',3);
    label4->setText(txt);
    poly2->Modified();
    qVTK2->GetRenderWindow()->Render();
}

double_t RoughnessWidget::calculateAsphericity(vtkSmartPointer<vtkPolyData> pd){
    auto N = pd->GetNumberOfPoints();
    VectorXd R(N);
    Vector3d p;
    double_t asph = 0;
    for(auto i=0; i < N; ++i){
        pd->GetPoint(i,&p(0));
        R(i) = p.norm();
    }
    asph = ((R.array()-R.mean()).square()).mean()/(R.mean()*R.mean());
    return asph;
}

double_t RoughnessWidget::calculateRMSAd(vtkSmartPointer<vtkPolyData> p){
    VectorXd Ad( p->GetNumberOfPoints() );
    Ad = 2*M_PI*VectorXd::Ones( p->GetNumberOfPoints() );
    auto pts = vtkSmartPointer<vtkIdList>::New();
    auto polys = p->GetPolys();
    polys->InitTraversal();
    while( polys->GetNextCell(pts) ){
        Vector3d p1, p2, p3, e1, e2, e3;
        p->GetPoint( pts->GetId(0), &p1(0) );
        p->GetPoint( pts->GetId(1), &p2(0) );
        p->GetPoint( pts->GetId(2), &p3(0) );
        e1 = (p3 - p2).normalized();
        e2 = (p1 - p3).normalized();
        e3 = (p2 - p1).normalized();
        Ad( pts->GetId(0) ) -= std::acos( e3.dot(-e2) );
        Ad( pts->GetId(1) ) -= std::acos( e1.dot(-e3) );
        Ad( pts->GetId(2) ) -= std::acos( e2.dot(-e1) );
    }
    double_t RMSad = std::sqrt( Ad.array().square().mean() );
    return RMSad;
}

void RoughnessWidget::on_slider1_valueChanged(double position){
    double_t d = 1 + 2.0*position;
    for(auto i=0; i < poly1->GetNumberOfPoints(); ++i){
        double p[3] = {0,0,0};
        poly1_0->GetPoint(i,p);
        p[1] = d*p[1];
        poly1->GetPoints()->SetPoint(i,p);
    }
    double_t val = calculateAsphericity(poly1);
    QString txt = QString("Asphericity: %1").arg(val,0,'g',3);
    label1->setText(txt);
    val = calculateRMSAd(poly1);
    txt = QString("RMS Angle Deficit: %1").arg(val,0,'g',3);
    label3->setText(txt);
    poly1->Modified();
    qVTK1->GetRenderWindow()->Render();
}
