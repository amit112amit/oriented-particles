#include "HelperFunctions.h"

namespace OPS{

void delaunay3DSurf(vtkSmartPointer<vtkPolyData> poly, std::string fi){
    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    auto pts = vtkSmartPointer<vtkPoints>::New();
    auto pts2 = vtkSmartPointer<vtkPoints>::New();
    auto unitSphere = vtkSmartPointer<vtkPolyData>::New();
    auto dssf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    auto idf = vtkSmartPointer<vtkIdFilter>::New();
    auto d3D = vtkSmartPointer<vtkDelaunay3D>::New();
    auto finalCells = vtkSmartPointer<vtkCellArray>::New();
    auto pointIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPolyData> final;
    vtkSmartPointer<vtkCellArray> interim;
    vtkSmartPointer<vtkIdTypeArray> origIds;

    // Read the VTK file
    reader->SetFileName(fi.c_str());
    reader->Update();
    size_t N = reader->GetOutput()->GetNumberOfPoints();
    pts = reader->GetOutput()->GetPoints();
    Eigen::Map<Eigen::Matrix3Xd> origP(
                            (double_t*)pts->GetData()->GetVoidPointer(0),3,N);
    Eigen::Matrix3Xd unitP(3,N);
    for(int i=0; i < N; i++){
        unitP.col(i) = origP.col(i).normalized();
        //We need to do this hack to get vtkDelaunay3D to behave well
        unitP(0,i) = std::round( unitP(0,i)*100 )/100;
        unitP(1,i) = std::round( unitP(1,i)*100 )/100;
        unitP(2,i) = std::round( unitP(2,i)*100 )/100;
    }
    void * voidPtr = (void*) unitP.data();
    auto unitPData = vtkSmartPointer<vtkDoubleArray>::New();
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
    return;
}

//! Calculate average edge length for a polydata
double_t getPointCloudAvgEdgeLen(std::string file){
    // Get a polydata mesh
    auto pd = vtkSmartPointer<vtkPolyData>::New();
    delaunay3DSurf(pd,file);
    // Now calculate average edge length
    double_t len = 0;
    auto extract = vtkSmartPointer<vtkExtractEdges>::New();
    extract->SetInputData(pd);
    extract->Update();
    auto edgePoly = extract->GetOutput();
    auto ptIds = vtkSmartPointer<vtkIdList>::New();
    auto lines = edgePoly->GetLines();
    lines->InitTraversal();
    Eigen::Vector3d p1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d p2 = Eigen::Vector3d::Zero();
    while( lines->GetNextCell(ptIds) ){
        edgePoly->GetPoint( ptIds->GetId(0), &p1(0) );
        edgePoly->GetPoint( ptIds->GetId(1), &p2(0) );
        len += (p1 - p2).norm();
    }
    len /= edgePoly->GetNumberOfLines();
    return len;
}

}
