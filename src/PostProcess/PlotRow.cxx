#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef Delaunay::Face_circulator Face_circulator;
typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Point Point;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix3Xd Matrix3Xd;
typedef Eigen::Map<Matrix3Xd> Map3Xd;
typedef Eigen::Quaterniond Quaterniond;
typedef Eigen::AngleAxisd AngleAxisd;

int main(int argc, char *argv[])
{

  // We expect four command line arguments:
  // 1. Number of particles in the simulation
  // 2. Name of data file from which to read rows
  // 3. Row number to plot -- 0 indexed row number
  // 4. Output VTK file name
  if (argc < 4)
  {
    std::cout << "Usage: ./plotRow <number_of_particles> <data_file> "
              << "<row_num> <output_vtk_file>" << std::endl;
    return 0;
  }

  // The number of particles in the simulation
  size_t N = std::stoi(std::string(argv[1]));

  // The row number to be written to vtk file
  size_t row = std::stoi(std::string(argv[3]));

  std::ifstream inputfile(argv[2]);
  if (!inputfile.is_open())
  {
    std::cout << "Unable to open file " << argv[2] << std::endl;
    return 0;
  }

  //**********************************************************************//
  // Now read the data file line by line
  //

  // The first line contains gamma and beta information
  std::string line, ignore;
  double_t gamma = 0, beta = 0;
  std::getline(inputfile, line);
  std::istringstream header(line);
  header >> ignore >> gamma >> ignore >> beta;

  // Skip to the required row
  auto rowCount = 0;
  while (rowCount < row)
  {
    std::getline(inputfile, line);
    rowCount++;
    continue;
  }

  std::getline(inputfile, line);
  std::cout << "Processing row " << rowCount << std::endl;
  std::istringstream rowStream(line);

  // We will extract 3*N doubles representing particle positions
  // from the stream
  std::string value;
  size_t valCount = 0;
  Matrix3Xd positions(3, N);
  while (valCount < 3 * N && rowStream.good())
  {
    std::getline(rowStream, value, ',');
    positions(valCount % 3, valCount / 3) = std::stod(value);
    valCount++;
  }

  // Read the rotation vectors
  Matrix3Xd rotVecs(3, N);
  valCount = 0;
  while (valCount < 3 * N && rowStream.good())
  {
    std::getline(rowStream, value, ',');
    rotVecs(valCount % 3, valCount / 3) = std::stod(value);
    valCount++;
  }

  // Calculate point normals using the rotation vectors
  Matrix3Xd normals(3, N);
  Quaterniond zaxis(0.0, 0.0, 0.0, 1.0);
  for (auto i = 0; i < N; ++i)
  {
    normals.col(i) = (Quaterniond(AngleAxisd(rotVecs.col(i).norm(),
                                             rotVecs.col(i).normalized())) *
                      zaxis *
                      (Quaterniond(AngleAxisd(rotVecs.col(i).norm(),
                                              rotVecs.col(i).normalized()))
                           .conjugate()))
                         .vec();
  }

  // Stereo graphic projection
  Matrix3Xd points(3, N);

  // Project points to unit sphere
  points = positions.colwise().normalized();

  // Reset the center of the sphere to origin by translating
  Vector3d center = points.rowwise().mean();
  points = points.colwise() - center;

  // Rotate all points so that the point in 0th column is along z-axis
  Vector3d c = points.col(0);
  double_t cos_t = c(2);
  double_t sin_t = std::sqrt(1 - cos_t * cos_t);
  Vector3d axis;
  axis << c(1), -c(0), 0.;
  Matrix3d rotMat, axis_cross, outer;
  axis_cross << 0., -axis(2), axis(1), axis(2), 0., -axis(0), -axis(1), axis(0),
      0.;

  outer.noalias() = axis * axis.transpose();

  rotMat =
      cos_t * Matrix3d::Identity() + sin_t * axis_cross + (1 - cos_t) * outer;
  Matrix3Xd rPts(3, N);
  rPts = rotMat * points; // The points on a sphere rotated

  // Calculate the stereographic projections
  Vector3d p0;
  Map3Xd l0(&(rPts(0, 1)), 3, N - 1);
  Matrix3Xd l(3, N - 1), proj(3, N - 1);
  p0 << 0, 0, -1;
  c = rPts.col(0);
  l = (l0.colwise() - c).colwise().normalized();
  for (auto j = 0; j < N - 1; ++j)
  {
    proj.col(j) = ((p0(2) - l0(2, j)) / l(2, j)) * l.col(j) + l0.col(j);
  }

  // Insert the projected points in a CGAL vertex_with_info vector
  std::vector<std::pair<Point, unsigned>> verts;
  for (auto j = 0; j < N - 1; ++j)
  {
    verts.push_back(std::make_pair(Point(proj(0, j), proj(1, j)), j + 1));
  }

  // Triangulate
  Delaunay dt;
  dt.insert(verts.begin(), verts.end());

  // Write the finite faces of the triangulation to a VTK file
  auto triangles = vtkSmartPointer<vtkCellArray>::New();
  for (auto ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end();
       ++ffi)
  {
    triangles->InsertNextCell(3);
    for (auto j = 2; j >= 0; --j)
      triangles->InsertCellPoint(ffi->vertex(j)->info());
  }

  // Iterate over infinite faces
  Face_circulator fc = dt.incident_faces(dt.infinite_vertex()), done(fc);
  if (fc != 0)
  {
    do
    {
      triangles->InsertNextCell(3);
      for (auto j = 2; j >= 0; --j)
      {
        auto vh = fc->vertex(j);
        auto id = dt.is_infinite(vh) ? 0 : vh->info();
        triangles->InsertCellPoint(id);
      }
    } while (++fc != done);
  }

  // Write to vtk file
  auto ptsArr = vtkSmartPointer<vtkDoubleArray>::New();
  ptsArr->SetVoidArray((void *)positions.data(), 3 * N, 1);
  ptsArr->SetNumberOfComponents(3);
  auto pts = vtkSmartPointer<vtkPoints>::New();
  pts->SetData(ptsArr);
  auto poly = vtkSmartPointer<vtkPolyData>::New();
  poly->SetPoints(pts);
  poly->SetPolys(triangles);
  auto normArr = vtkSmartPointer<vtkDoubleArray>::New();
  normArr->SetVoidArray((void *)normals.data(), 3 * N, 1);
  normArr->SetNumberOfComponents(3);
  poly->GetPointData()->SetNormals(normArr);
  auto wr = vtkSmartPointer<vtkPolyDataWriter>::New();
  wr->SetFileName(argv[4]);
  wr->SetInputData(poly);
  wr->Write();

  inputfile.close();
  return 0;
}
