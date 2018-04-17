/* Reference:
 *       Szeliski, Richard and Tonnesen, D. (1992). Surface Modeling with
 *       Oriented Particle Systems.	Siggraph ’92, 26(2), 160.
 *       https://doi.org/10.1017/CBO9781107415324.004
 */
#if !defined(__OPSBODY_H__)
#define __OPSBODY_H__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkOctreePointLocator.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractEdges.h>
#include "Body.h"
#include "HelperFunctions.h"

namespace OPS{
    // ************************* OPSBody Class *********************************//
    /*!
     * \brief The OPSBody class
     * Computes energy and forces on a Oriented Particle System
     */
    class OPSBody : public Body{
	public:
	    typedef Eigen::Vector3d Vector3d;
	    typedef Eigen::Map<Vector3d> MapV3d;
	    typedef Eigen::Ref<const Vector3d> RefCV3d;
	    typedef Eigen::VectorXd VectorXd;
	    typedef Eigen::Matrix3d Matrix3d;
	    typedef Eigen::Ref<Matrix3d> RefM3d;
	    typedef Eigen::Matrix3Xd Matrix3Xd;
	    typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
	    typedef Eigen::Map<Matrix3Xd> MapM3Xd;
	    typedef Eigen::Matrix<double_t, 3, 4> Matrix3x4d;
	    typedef Eigen::Matrix<double_t, 4, 3> Matrix4x3d;
	    typedef Eigen::Quaterniond Quaterniond;
	    typedef Eigen::AngleAxisd AngleAxisd;

	    OPSBody(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot, RefM3Xd posGrad,
		    RefM3Xd rotGrad);
	    void applyKabschAlgorithm();
	    virtual void compute();
	    void computeNormals();
	    void diffNormalRotVec();
	    double_t determineSearchRadius();
	    double_t getAsphericity();
	    double_t getArea();
	    double_t getAverageEdgeLength();
	    double_t getAverageNumberOfNeighbors(){return _avgNumNeighbors;}
	    double_t getAverageRadius();
	    double_t getCircularityEnergy(){return _circEn;}
	    double_t getMeanSquaredDisplacement();
	    std::vector<double_t> getMSD();
	    double_t getMorseEnergy(){return _morseEn;}
	    double_t getNormalityEnergy(){return _normalEn;}
	    vtkSmartPointer<vtkPolyData> getPolyData(){return _polyData;}
	    double_t getRMSAngleDeficit();
	    double_t getTotalEnergy();
	    double_t getVolume();
	    static void initialRotationVector(RefM3Xd pos, RefM3Xd rotVec);
	    void printVTKFile(const std::string name);
	    void resetToInitialPositions();
	    void saveInitialPosition(){ _initialPositions = _positions;}
	    void setMorseDistance(double_t r){_re = r;}
	    void setMorseWellWidth(double_t a){_a = a;}
	    void setSearchRadius(double_t s){_searchRadius = s;}
	    void setFVK(double_t g){_gamma = g;}
	    void sphericalDelaunay();
	    virtual void updateNeighbors();
	    void updateDataForKabsch();
	    void updatePolyData();

	protected:
	    bool _updateArea = true, _updateRadius = true, _updateVolume = true;
	    double_t &_f, _area = 0, _radius = 0, _volume = 0, _morseEn, _normalEn,
		     _circEn, _msd_tgt = 0, _avgNumNeighbors = 0, _re = 1.0,
		     _a = 4.62098120373, _gamma = 0.1, _searchRadius = 1.4,
		     _msd;
	    int _N;
	    MapM3Xd _positions, _posGradient, _rotGradient, _rotationVectors;
	    Matrix3Xd _prevX, _normals, _initialPositions;
	    std::vector< Matrix3d > _diffNormalRV; /*!< Derivatives of normals with rotation vectors*/
	    std::vector< vtkSmartPointer<vtkIdList> > _neighbors;
	    std::vector< vtkIdType > _initialNearestNeighbor;
        vtkSmartPointer< vtkPolyData > _polyData;
	    vtkSmartPointer< vtkOctreePointLocator > _octree;
	    void updateRotationVectors();
    };
    //***************************************************************************//
}
#endif //__OPSBODY_H__
