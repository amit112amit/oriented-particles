#if !defined(__OPSMODEL_H__)
#define __OPSMODEL_H__

#include <stdio.h>
#include <set>
#include <array>
#include <string.h>
#include <math.h>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vtkSmartPointer.h>
#include <vtkOctreePointLocator.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkIdFilter.h>
#include <vtkPointData.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "HelperFunctions.h"
#include "Body.h"

namespace OPS{

typedef Eigen::Vector3d Vector3d;

struct vertex_first_ring{
    size_t vertex_id;
    std::vector< std::pair<unsigned,unsigned> > faces;
    std::vector< unsigned > edges;
};

class OPSModel: public Body{
public:
        typedef Eigen::Map<Vector3d> MapV3d;
        typedef Eigen::Ref<const Vector3d> RefCV3d;
        typedef Eigen::VectorXd VectorXd;
        typedef Eigen::Ref<VectorXd> RefVXd;
        typedef Eigen::Matrix3d Matrix3d;
        typedef Eigen::Ref<Matrix3d> RefM3d;
        typedef Eigen::Matrix3Xd Matrix3Xd;
        typedef Eigen::Ref<Matrix3Xd> RefM3Xd;
        typedef Eigen::Map<Matrix3Xd> MapM3Xd;
        typedef Eigen::Matrix<double_t, 3, 4> Matrix3x4d;
        typedef Eigen::Matrix<double_t, 4, 3> Matrix4x3d;
        typedef Eigen::Quaterniond Quaterniond;
        typedef Eigen::AngleAxisd AngleAxisd;
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_2 Point;
        typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned,K> Vb;
        typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
        typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;
        typedef Delaunay::Face_circulator Face_circulator;

        OPSModel(size_t n, double_t &f, RefVXd x, RefVXd g, RefVXd pX);
        void applyKabschAlgorithm();
        virtual void compute();
        void computeNormals();
        void diffNormalRotVec();
        double_t getAsphericity();
        double_t getArea();
        double_t getAverageEdgeLength();
        double_t getAverageRadius();
        inline double_t getCircularityEnergy(){return _circEn;}
        inline void getInitialPositions(Matrix3Xd &v){ v = _initialPositions;}
        double_t getMSD();
        inline std::vector<size_t> getInitialNeighbors(){
                return _initialNearestNeighbor;
        }
        inline double_t getMorseEnergy(){return _morseEn;}
        inline double_t getNormalityEnergy(){return _normalEn;}
        double_t getRMSAngleDeficit();
        inline double_t getTotalEnergy(){ return (_morseEn + _normalEn + _circEn);}
        double_t getVolume();
        static void getRotationVectors(RefM3Xd pos, RefM3Xd rotVec);
        void printVTKFile(const std::string name);
        void resetToInitialPositions();
        inline void saveInitialPosition(){ _initialPositions = _positions;}
        inline void setInitialNeighbors(std::vector<size_t> &x){_initialNearestNeighbor = x;}
        inline void setInitialPositions(Matrix3Xd p){_initialPositions = p;}
        inline void setMorseDistance(double_t r){_re = r;}
        inline void setMorseWellWidth(double_t a){_a = a;}
        inline void setFVK(double_t g){_gamma = g;}
        void stereoDelaunay();
        void updateTriangles();
        inline bool constraintSatisfied(){
                return (std::abs(_value - _constrainedValue) < _tolerance);
        }
        inline void setConstraint(double_t V){_constrainedValue = V;}
        inline void setLagrangeCoeff(double_t L){_Lambda_i = L;}
        inline void setPenaltyCoeff(double_t K){_K_i = K;}
        inline void setTolerance(double_t t){_tolerance = t;}
        inline void uzawaUpdate(){
                _Lambda_i = _Lambda_i - _K_i*(_value - _constrainedValue);
                _K_i *= 10;
        }
protected:
        bool _updateArea = true, _updateRadius = true, _updateVolume = true;
        double_t _area = 0, _radius = 0, _volume = 0, _morseEn, _normalEn,
        _circEn, _re = 1.0, _a = 4.62098120373, _gamma = 0.1, _msd;
        int _N;
        MapM3Xd _positions, _posGradient, _rotGradient, _rotationVectors, _prevX;
        Matrix3Xd _normals, _initialPositions;
        std::vector< Matrix3d, Eigen::aligned_allocator<Matrix3d> > _diffNormalRV;
        std::vector<size_t> _initialNearestNeighbor;
        double_t &_f;
        double_t _constrainedValue = 0.0;
        double_t _K_i = 1000.0;
        double_t _Lambda_i = 1.0;
        double_t _tolerance = 1e-10;
        double_t _value = 0.0;
        Delaunay _dt;
        std::vector<vertex_first_ring> _first_ring;
};

}

#endif //__OPSMODEL_H__
