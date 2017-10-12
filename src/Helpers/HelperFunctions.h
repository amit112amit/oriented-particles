#if !defined(__HELPERFUNCTIONS_H__)
#define __HELPERFUNCTIONS_H__

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

//! Rotates a point cloud to get rid of rigid body motions
Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
                                      Eigen::Ref<Eigen::Matrix3Xd> out);

//! Generates a mesh for a topologically spherical point cloud and writes
//! it to file
void delaunay3DSurf(vtkSmartPointer<vtkPoints> pts, std::string fileName) ;

#endif // __HELPERFUNCTIONS_H__