
#ifndef SRC_MESH_ELEMENTS_POINT_DOF_H_
#define SRC_MESH_ELEMENTS_POINT_DOF_H_

#include "pointelement.h"

#define DOFGPCount 0
#define DOFVTKCode -4

namespace espreso {

class DOF: public PointElement
{

public:
	DOF(eslocal index): PointElement(index) {};
	Element* copy() const { return new DOF(*this); }

	eslocal vtkCode() const { return DOFVTKCode; }
	size_t gaussePoints() const { return DOFGPCount; }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }

	const std::vector<Property>& elementDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& faceDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& edgeDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& pointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& midPointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
};

}



#endif /* SRC_MESH_ELEMENTS_POINT_DOF_H_ */
