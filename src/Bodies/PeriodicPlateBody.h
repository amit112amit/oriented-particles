#if !defined(__MORSE2DBODY_H__)
#define __MORSE2DBODY_H__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include "Body.h"
#include "HelperFunctions.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkKdTreePointLocator.h>

namespace OPS{
    // *********************** OPSBody Class *******************************//
    /*!
     * \brief The Morse2DBody class
     * Computes energy and forces on a 2D Hexagonal lattice interacting with Morse
     * potentials
     */
    class PeriodicPlateBody : public Body{
	public:
	    typedef Eigen::Vector3d Vector3d;
	    typedef Eigen::Matrix3d Matrix3d;
	    typedef Eigen::Matrix3Xd Matrix3Xd;
	    typedef Eigen::Map<Matrix3Xd> Map3Xd;
	    PeriodicPlateBody(size_t N, double_t &f, Matrix3Xd &p,
		    Matrix3Xd &g, Matrix3Xd &r, Matrix3Xd &rg, double_t L,
		    double_t W);
	    void applyKabschAlgorithm();
	    virtual void compute();
	    void computeNormals();
	    void diffNormalRotVec();
	    double_t getCircularityEnergy(){return _circEn;}
	    double_t getMeanSquaredDisplacement();
	    double_t getMorseEnergy(){return _morseEn;}
	    double_t getNormalityEnergy(){return _normEn;}
	    static void initialRotationVector(Matrix3Xd &pos,
		    Matrix3Xd &rotVec);
	    void printVTKFile(const std::string name);
	    void setCircularityCoeff(double_t v){ _circCoeff = v;}
	    void setFVK(double_t v){_gamma = v;}
	    void setMorseDistance(double_t r){_re = r;}
	    void setMorseWellWidth(double_t a){_a = a;}
	    void setSearchRadius(double_t s){_searchRadius = s;}
	    void updateDataForKabsch();
	    void updatePeriodicKdTree();
	    void updateNeighbors();

	protected:
	    double_t &_f, _morseEn, _normEn, _circEn, _msd = 0, _De = 1.0,
		     _re = 1.0, _a = 4.62098120373, _searchRadius = 1.4,
		     _L = 10.0, _W = 10.0, _gamma = 0.1, _circCoeff = 1.0;
	    size_t _N;
	    Map3Xd _positions, _posGradient, _rotGradient, _rotationVectors;
	    Matrix3Xd  _initialPositions, _normals, _prevX;
	    std::vector< Matrix3d > _diffNormalRV;
	    std::vector< vtkSmartPointer<vtkIdList> > _neighbors;
	    std::vector< vtkIdType > _initialNearestNeighbor;
	    std::map< vtkIdType, vtkIdType > periodicMap;
	    vtkSmartPointer< vtkPolyData > _polyData;
	    vtkSmartPointer< vtkKdTreePointLocator > _kdtree;
	    void updateRotationVectors();

    };
    //***********************************************************************//
}
#endif //__MORSE2DBODY_H__
