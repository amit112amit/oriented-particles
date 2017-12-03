#if !defined(__ALCONSTRAINT_H__)
#define __ALCONSTRAINT_H__

#include <Eigen/Dense>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include "Body.h"

//! Parent class for all Augmented Lagrangian constraints
class ALConstraint: public Body{
	public:
		typedef Eigen::Matrix3Xd Matrix3Xd;
		typedef Eigen::Map<Matrix3Xd> MapM3Xd;
		typedef Eigen::Ref<Matrix3Xd> RefM3Xd;

		ALConstraint(size_t N, double_t &f, RefM3Xd x, RefM3Xd g);
		virtual void compute() = 0;
		virtual bool constraintSatisfied(){
			return (std::abs(_value - _constrainedValue) < _tolerance);
		}
		virtual void printCompletion(){
			std::cout<< " Constraint absolute value = "
				<< std::abs(_value - _constrainedValue) << std::endl;
		}
		void setConstraint(double_t V){_constrainedValue = V;}
		void setLagrangeCoeff(double_t L){_Lambda_i = L;}
		void setPenaltyCoeff(double_t K){_K_i = K;}
		void setTolerance(double_t t){_tolerance = t;}
		void uzawaUpdate();

	protected:
		double_t _constrainedValue = 0.0;
		double_t &_f;
		double_t _K_i = 1000.0;
		double_t _Lambda_i = 1.0;
		double_t _tolerance = 1e-10;
		double_t _value = 0.0;
		MapM3Xd _xPos;
		MapM3Xd _xGrad;
		size_t _N;
};

//! Implements an average area Augmented Lagrangian constraint
class AvgAreaConstraint: public ALConstraint{
	public:
		AvgAreaConstraint(size_t N, double_t &f, RefM3Xd x,
				RefM3Xd g):ALConstraint(N,f,x,g){}
		void compute();
};

//! Implements an average volume Augmented Lagrangian constraint
class AvgVolConstraint: public ALConstraint{
	public:
		AvgVolConstraint(size_t N, double_t &f, RefM3Xd x,
				RefM3Xd g):ALConstraint(N,f,x,g){}
		void compute();
};

//! Implements an exact area Augmented Lagrangian constraint
class ExactAreaConstraint: public ALConstraint{
	public:
		ExactAreaConstraint(size_t N, double_t &f, RefM3Xd x,
				RefM3Xd g, vtkSmartPointer<vtkPolyData> p):
			ALConstraint(N,f,x,g){
				_poly = p;
			}
		void compute();
	private:
		vtkSmartPointer<vtkPolyData> _poly;
};

//! Implements an exact volume Augmented Lagrangian constraint
class ExactVolConstraint: public ALConstraint{
	public:
		ExactVolConstraint(size_t N, double_t &f, RefM3Xd x,
				RefM3Xd g, vtkSmartPointer<vtkPolyData> p):
			ALConstraint(N,f,x,g){
				_poly = p;
			}
		void compute();
	private:
		vtkSmartPointer<vtkPolyData> _poly;
};

//! Implements both an exact area and an exact volume constraint
class ExactAreaVolConstraint: public ALConstraint{
	public:
		ExactAreaVolConstraint(size_t N, double_t &f, RefM3Xd x,
				RefM3Xd g, vtkSmartPointer<vtkPolyData> p):
			ALConstraint(N,f,x,g){
				_poly = p;
				_area = 0.0;
				_volume = 0.0;
				_areaConstrained = 0.0;
				_volConstrained = 0.0;
			}    
		bool constraintSatisfied(){
			return (std::abs(_area - _areaConstrained) < _tolerance) &&
				(std::abs(_volume - _volConstrained) < _tolerance);
		}
		void printCompletion(){
			std::cout<<"Area constraint = "<<std::abs(_area - _areaConstrained)
				<<" Volume constraint = "<<std::abs(_volume - _volConstrained)
				<< std::endl;
		}
		using ALConstraint::setConstraint;
		void setConstraint(double_t A, double_t V){
			_areaConstrained = A;
			_volConstrained = V;
		}
		void compute();
	private:
		double_t _area;
		double_t _volume;
		double_t _volConstrained;
		double_t _areaConstrained;
		vtkSmartPointer<vtkPolyData> _poly;
};
#endif //__ALCONSTRAINT_H__
