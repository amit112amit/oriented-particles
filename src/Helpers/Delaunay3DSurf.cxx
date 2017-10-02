#include "HelperFunctions.h"

void delaunay3DSurf(vtkSmartPointer<vtkPoints> pts, std::string fileName ) {
    vtkSmartPointer<vtkPoints> pts2 =
            vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> unitSphere =
            vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> final;
    vtkSmartPointer<vtkDataSetSurfaceFilter> dssf =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    vtkSmartPointer<vtkIdFilter> idf =
            vtkSmartPointer<vtkIdFilter>::New();
    vtkSmartPointer<vtkDelaunay3D> d3D =
            vtkSmartPointer<vtkDelaunay3D>::New();
    vtkSmartPointer<vtkCellArray> finalCells =
            vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> interim;
    vtkSmartPointer<vtkIdTypeArray> origIds;
    vtkSmartPointer<vtkIdList> pointIds =
            vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPolyData> poly =
            vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyDataWriter> w =
            vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkKdTree> kdTree =
            vtkSmartPointer<vtkKdTree>::New();

    kdTree->BuildLocatorFromPoints( pts );

    size_t N = pts->GetNumberOfPoints();
    Eigen::Map<Eigen::Matrix3Xd> origP(
                (double_t*)pts->GetData()->GetVoidPointer(0),
                3,N);
    Eigen::Matrix3Xd unitP(3,N);
    for(int i=0; i < N; i++){
        unitP.col(i) = origP.col(i).normalized();
        //We need to do this hack to get vtkDelaunay3D to behave well
        unitP(0,i) = std::round( unitP(0,i)*100 )/100;
        unitP(1,i) = std::round( unitP(1,i)*100 )/100;
        unitP(2,i) = std::round( unitP(2,i)*100 )/100;
    }
    void * voidPtr = (void*) unitP.data();
    vtkSmartPointer<vtkDoubleArray> unitPData =
            vtkSmartPointer<vtkDoubleArray>::New();
    unitPData->SetVoidArray(voidPtr,3*N,1);
    unitPData->SetNumberOfComponents(3);
    pts2->SetData(unitPData);

    unitSphere->SetPoints(pts2);

    idf->SetIdsArrayName("PointIds");
    idf->PointIdsOn();
    idf->SetInputData(unitSphere);

    d3D->SetInputConnection(idf->GetOutputPort());
    dssf->SetInputConnection(d3D->GetOutputPort());
    dssf->Update();
    final = dssf->GetOutput();
    interim = final->GetPolys();
    interim->InitTraversal();
    origIds = vtkIdTypeArray::SafeDownCast(
                final->GetPointData()->GetArray("PointIds"));
    while(interim->GetNextCell(pointIds)){
        int numIds = pointIds->GetNumberOfIds();
        finalCells->InsertNextCell(numIds);
        for(int j=0; j < numIds; j++ ){
            int id = (int)origIds->GetTuple1( pointIds->GetId(j) );
            finalCells->InsertCellPoint(id);
        }
    }

    poly->SetPoints(pts);
    poly->SetPolys(finalCells);
    w->SetInputData(poly);
    w->SetFileName(fileName.c_str());
    w->Write();
}
