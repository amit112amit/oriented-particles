#include "OPSMesh.h"

namespace OPS
{
//! Constructor for OPSMesh
OPSMesh::OPSMesh(size_t n, double_t &f, double_t R0, RefM3Xd pos, RefM3Xd rot,
                 RefM3Xd pG, RefM3Xd rG, RefM3Xd pX)
    : OPSBody(n, f, R0, pos, rot, pG, rG, pX)
{
  _edges = vtkSmartPointer<vtkCellArray>::New();
  _edgePoly = vtkSmartPointer<vtkPolyData>::New();
  _numBonds = (int)((12 * 5 + (n - 12) * 6) / 2);
}

//! Extract the edges of the polydata
//! assuming that the _polyData has been updated
void OPSMesh::updateNeighbors()
{
  // Create the VTK objects
  auto extract = vtkSmartPointer<vtkExtractEdges>::New();
  auto idf = vtkSmartPointer<vtkIdFilter>::New();
  auto newEdges = vtkSmartPointer<vtkCellArray>::New();
  auto origIds = vtkSmartPointer<vtkIdTypeArray>::New();
  auto verts = vtkSmartPointer<vtkIdList>::New();

  // Extract the edges
  idf->SetPointIdsArrayName("PointIds");
  idf->PointIdsOn();
  idf->SetInputData(_polyData);
  extract->SetInputConnection(idf->GetOutputPort());
  extract->Update();
  newEdges = extract->GetOutput()->GetLines();
  origIds = vtkIdTypeArray::SafeDownCast(
      extract->GetOutput()->GetPointData()->GetArray("PointIds"));

  // Convert edges in terms of the original point ids
  _edges->Reset();
  newEdges->InitTraversal();
  while (newEdges->GetNextCell(verts))
  {
    _edges->InsertNextCell(2);
    vtkIdType id = (vtkIdType)origIds->GetTuple1(verts->GetId(0));
    _edges->InsertCellPoint(id);
    id = (vtkIdType)origIds->GetTuple1(verts->GetId(1));
    _edges->InsertCellPoint(id);
  }
}

//! Compute the OPSBody energy
void OPSMesh::compute()
{
  auto pts = vtkSmartPointer<vtkIdList>::New();

  auto factor = 0.666666666666667 * _a * _a * _R0 * _R0 / _gamma;

  // Initialize energies and forces to be zero
  _morseEn = 0.0;
  _normalEn = 0.0;
  _circEn = 0.0;

  computeNormals();
  diffNormalRotVec();

  _edges->InitTraversal();
  while (_edges->GetNextCell(pts))
  {

    double_t r, n_dot_rn, exp_1, exp_2, morseEn, Phi_n, Phi_c;
    Matrix3d M, N;
    Vector3d vi, p, vj, q, m, n, rij, rn, dMdr, dPhi_nVi, dPhi_nVj, dPhi_cVi,
        dPhi_cVj, dCdr, Dxi, Dvi, Dvj;
    Vector3d xi = Vector3d::Zero();
    Vector3d xj = Vector3d::Zero();
    vtkIdType i, j;

    i = pts->GetId(0);
    j = pts->GetId(1);

    xi = _positions.col(i);
    xj = _positions.col(j);
    vi = _rotationVectors.col(i);
    vj = _rotationVectors.col(j);
    p = _normals.col(i);
    q = _normals.col(j);
    M = _diffNormalRV[i];
    N = _diffNormalRV[j];

    rij = xj - xi;
    rn = rij.normalized();
    r = rij.norm();
    m = p - q;
    n = p + q;
    n_dot_rn = n.dot(rn);

    // Evaluate morse derivatives
    exp_1 = exp(-_a * (r - _re));
    exp_2 = exp_1 * exp_1;
    morseEn = exp_2 - 2 * exp_1;
    dMdr = 2 * _a * (exp_1 - exp_2) * rn;

    // Evaluate co-normality derivatives
    Phi_n = m.squaredNorm();
    dPhi_nVi = 2 * M * m;
    dPhi_nVj = -2 * N * m;

    // Evaluate co-circularity derivatives
    Phi_c = n_dot_rn * n_dot_rn;
    dCdr = (2 * n_dot_rn / r) * (n - n_dot_rn * rn);
    dPhi_cVi = (2 * n_dot_rn) * M * rn;
    dPhi_cVj = (2 * n_dot_rn) * N * rn;

    // Calculate the total derivatives of energy wrt xi, vi and vj
    Dxi = -(dMdr + dCdr * factor);
    Dvi = (dPhi_nVi + dPhi_cVi) * factor;
    Dvj = (dPhi_nVj + dPhi_cVj) * factor;

    // Update the energies
    _morseEn += morseEn;
    _normalEn += Phi_n * factor;
    _circEn += Phi_c * factor;
    _f += morseEn + (Phi_n + Phi_c) * factor;

    // Update the derivatives
    _posGradient.col(i) += Dxi;
    _posGradient.col(j) -= Dxi; /*!< Neighbor feels opposite force */
    _rotGradient.col(i) += Dvi;
    _rotGradient.col(j) += Dvj;
  }
  return;
}

//! Function to return the difference between normals of all points connected
//! with a bond
void OPSMesh::getDiffNormals(Matrix3Xd &in)
{
  in = Matrix3Xd::Zero(3, _numBonds);
  auto pts = vtkSmartPointer<vtkIdList>::New();
  _edges->InitTraversal();
  auto index = 0;
  while (_edges->GetNextCell(pts))
  {
    Vector3d p, q;
    vtkIdType i, j;

    i = pts->GetId(0);
    j = pts->GetId(1);

    p = _normals.col(i);
    q = _normals.col(j);
    in.col(index++) = p - q;
  }
}

} // namespace OPS
