#if !defined(__HELPERFUNCTIONS_H__)
#define __HELPERFUNCTIONS_H__

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <map>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

namespace OPS{

    //Structure declaration
    struct neighbors {
	vtkIdType _id;
	double _distance;
	neighbors(vtkIdType id, double dist) :_id(id), _distance(dist) {}
	bool operator<(const neighbors& e)const { return _distance < e._distance; }
    };

    //Unique cells
    struct uniqueTriangle {
	vtkIdType _id1, _id2, _id3;
	uniqueTriangle(vtkIdType id1, vtkIdType id2, 
		vtkIdType id3):_id1(id1), _id2(id2), _id3(id3) {}
	bool operator<(const uniqueTriangle& e)const { 
	    std::vector<vtkIdType> set1 = {_id1, _id2, _id3};
	    std::vector<vtkIdType> set2 = {e._id1, e._id2, e._id3};
	    std::sort(set1.begin(), set1.end());
	    std::sort(set2.begin(), set2.end());
	    if (set1[0] != set2[0]) 
		return (set1[0] < set2[0]);
	    else if (set1[1] != set2[1]) 
		return (set1[1] < set2[1]);
	    else 
		return (set1[2] < set2[2]);
	}
    };

    typedef std::map< std::string, std::string > InputParameters;

    //! Rotates a point cloud to get rid of rigid body motions
    Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
	    Eigen::Ref<Eigen::Matrix3Xd> out);

    //! Generates a mesh for a topologically spherical point cloud and writes
    //! it to file
    void delaunay3DSurf(vtkSmartPointer<vtkPoints> pts, std::string fileName) ;

    //! Reads a key-value input file and returns a map object
    InputParameters readKeyValueInput( std::string fileName );

}
#endif // __HELPERFUNCTIONS_H__
