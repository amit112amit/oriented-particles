#if !defined(__OPSMESH_H__)
#define __OPSMESH_H__

#include <vtkExtractEdges.h>
#include "OPSBody.h"

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
    void compute();
    void updateNeighbors();
private:
    vtkSmartPointer<vtkPolyData> _edgePoly;
    vtkSmartPointer<vtkCellArray> _edges;
};
//***************************************************************************//
#endif //__OPSMESH_H__
