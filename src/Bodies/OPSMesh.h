#if !defined(__OPSMESH_H__)
#define __OPSMESH_H__

#include <vtkExtractEdges.h>
#include "OPSBody.h"

namespace OPS{
// ************************* OPSMesh Class *********************************//
/*!
 * \brief The OPSMesh class
 * Computes energy and forces on a Oriented Particle System using the
 * triangulation
 */
class OPSMesh : public OPSBody{
	public:
		OPSMesh(size_t n, double_t &f, RefM3Xd pos, RefM3Xd rot, RefM3Xd posGrad,
				RefM3Xd rotGrad);
		void getDiffNormals(Matrix3Xd &in);
		void getNormals(Matrix3Xd &in);
		void compute();
		void updateNeighbors();
		void setSpontaneousCurvature(double_t);
	private:
		vtkSmartPointer<vtkPolyData> _edgePoly;
		vtkSmartPointer<vtkCellArray> _edges;
		bool _spontaneousCurvatureOn = false;
		double_t _spontaneousCurvature = 0.0;
		size_t _numBonds;
};
//***************************************************************************//
}
#endif //__OPSMESH_H__

