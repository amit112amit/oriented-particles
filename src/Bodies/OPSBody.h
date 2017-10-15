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
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include "Body.h"
#include "HelperFunctions.h"

// **************************** OPSParams *********************************//
/*!
 * \brief The OPSParams struct
 *
 * Encapsulates the various parameters needed by a BrownOPS class
 */
class OPSParams{
public:
    enum Parameter{D_eV, r_eV, aV, bV, searchRadiusV, gammaV};

    //! Initialize with default values
    OPSParams():D_e(1.0), r_e(1.0), a(6.9314718056),//Fracture strain = 10%
        b(1.0), gamma(1.0),
        searchRadius(1.4){}

    //! Initialize from double_ts
    OPSParams( double_t E, double_t r, double_t av, double_t bv, double_t sr,
               double_t s): D_e(E), r_e(r), a(av), b(bv), searchRadius(sr),
        gamma(s){}

    //! Accessor methods
    double getMorseEquilibriumEnergy(){return D_e;}
    double getMorseEquilibriumDistance(){return r_e;}
    double getMorseWellWidth(){return a;}
    double getKernelStdDev(){return b;}
    double getSearchRadius(){return searchRadius;}
    double getFVK(){return gamma;}

    //! Update a parameter value
    void updateParameter(Parameter p, double_t val);

private:
    double_t D_e, r_e, a; /*!< Morse potential parameters */
    double_t b;           /*!< Standard deviation of the Kernel Gaussian */
    double_t searchRadius;
    double_t gamma; /*!< Energy scale ratio */
};
// **************************************************************************//

// ************************* BrownOPS Class *********************************//
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
    typedef Eigen::Map<Matrix3Xd> Map3Xd;
    typedef Eigen::Matrix<double_t, 3, 4> Matrix3x4d;
    typedef Eigen::Matrix<double_t, 4, 3> Matrix4x3d;
    typedef Eigen::Quaterniond Quaterniond;
    typedef Eigen::AngleAxisd AngleAxisd;

    OPSBody(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot, RefM3Xd posGrad,
             RefM3Xd rotGrad, OPSParams &p);
    void applyKabschAlgorithm();
    void compute();    
    void computeNormals();    
    void diffNormalRotVec(const RefCV3d &vi, RefM3d diff);
    double_t getAsphericity();
    double_t getArea();
    double_t getAverageEdgeLength();
    double_t getAverageRadius();
    double_t getAverageNumberOfNeighbors();
    double_t getCircularityEnergy(){return _circEn;}
    double_t getMeanSquaredDisplacement();
    double_t getMorseEnergy(){return _morseEn;}
    double_t getNormalityEnergy(){return _normalEn;}
    vtkSmartPointer<vtkPolyData> getPolyData(){return _polyData;}
    double_t getTotalEnergy();
    double_t getVolume();
    static void initialRotationVector(RefM3Xd pos, RefM3Xd rotVec);
    void printVTKFile(const std::string name);
    void saveInitialPosition(){ _initialPositions = _positions;}           
    void updateNeighbors();
    void updateDataForKabsch();
    void updatePolyDataAndKdTree();    

private:
    bool _updateArea;
    bool _updateRadius;
    bool _updateVolume;
    double_t &_f;
    double_t _area;
    double_t _radius;
    double_t _volume;
    double_t _volConstrained;
    double_t _morseEn;
    double_t _normalEn;
    double_t _circEn;
    double_t _msd;
    int _numPartilces;
    Map3Xd _positions;
    Map3Xd _posGradient;
    Map3Xd _rotGradient;
    Map3Xd _rotationVectors;
    Matrix3Xd _prevX;
    Matrix3Xd _normals;    
    Matrix3Xd _initialPositions;
    OPSParams &_params;
    std::vector< vtkSmartPointer<vtkIdList> > _neighbors;
    std::vector< vtkIdType > _initialNearestNeighbor;
    vtkSmartPointer< vtkPolyData > _polyData;
    vtkSmartPointer< vtkKdTree > _kdTree;
    void updateRotationVectors();
};
//***************************************************************************//
#endif //__OPSBODY_H__
