#include "HelperFunctions.h"

namespace OPS
{

void delaunay3DSurf(vtkSmartPointer<vtkPolyData> poly, std::string fi)
{
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
      (double_t *)pts->GetData()->GetVoidPointer(0), 3, N);
  Eigen::Matrix3Xd unitP(3, N);
  for (int i = 0; i < N; i++)
  {
    unitP.col(i) = origP.col(i).normalized();
    // We need to do this hack to get vtkDelaunay3D to behave well
    unitP(0, i) = std::round(unitP(0, i) * 100) / 100;
    unitP(1, i) = std::round(unitP(1, i) * 100) / 100;
    unitP(2, i) = std::round(unitP(2, i) * 100) / 100;
  }
  void *voidPtr = (void *)unitP.data();
  auto unitPData = vtkSmartPointer<vtkDoubleArray>::New();
  unitPData->SetVoidArray(voidPtr, 3 * N, 1);
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
  origIds =
      vtkIdTypeArray::SafeDownCast(final->GetPointData()->GetArray("PointIds"));
  while (interim->GetNextCell(pointIds))
  {
    int numIds = pointIds->GetNumberOfIds();
    finalCells->InsertNextCell(numIds);
    for (int j = 0; j < numIds; j++)
    {
      int id = (int)origIds->GetTuple1(pointIds->GetId(j));
      finalCells->InsertCellPoint(id);
    }
  }

  poly->SetPoints(pts);
  poly->SetPolys(finalCells);
}

Delaunay delaunayStereo(Eigen::Ref<Eigen::Matrix3Xd> points)
{
  // Make stereographic projections
  Eigen::Matrix3Xd proj = stereographicProjection(points);

  // Insert the projected points in a CGAL vertex_with_info vector
  std::vector<std::pair<Point, unsigned>> verts;
  for (auto j = 0; j < proj.cols(); ++j)
    verts.push_back(std::make_pair(Point(proj(0, j), proj(1, j)), j + 1));

  // CGAL Delauany 2d triangulation
  Delaunay dt(verts.begin(), verts.end());
  return dt;
}

//! Calculate average edge length for a polydata
double_t getPointCloudAvgEdgeLen(std::string file)
{
  // Read the vtk file
  auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(file.c_str());
  reader->Update();
  auto N = reader->GetOutput()->GetNumberOfPoints();
  auto *arrPointer =
      (double_t *)(reader->GetOutput()->GetPoints()->GetData()->GetVoidPointer(
          0));
  Eigen::Map<Eigen::Matrix3Xd> positions(arrPointer, 3, N);

  // Get a triangulation mesh
  Delaunay dt = delaunayStereo(positions);

  // Now calculate average edge length
  std::vector<double_t> edgeLengths;
  for (auto fei = dt.finite_edges_begin(); fei != dt.finite_edges_end();
       ++fei)
  {
    unsigned edgeVert1, edgeVert2;
    auto fh = fei->first;
    switch (fei->second)
    {
    case 0:
      edgeVert1 = fh->vertex(1)->info();
      edgeVert2 = fh->vertex(2)->info();
      break;
    case 1:
      edgeVert1 = fh->vertex(2)->info();
      edgeVert2 = fh->vertex(0)->info();
      break;
    case 2:
      edgeVert1 = fh->vertex(0)->info();
      edgeVert2 = fh->vertex(1)->info();
      break;
    }
    edgeLengths.push_back(
        (positions.col(edgeVert1) - positions.col(edgeVert2)).norm());
  }
  auto avgEdgeLen =
      std::accumulate(edgeLengths.begin(), edgeLengths.end(), 0.0) /
      edgeLengths.size();

  return avgEdgeLen;
}

} // namespace OPS
