
#ifndef SRC_MESH_ELEMENTS_POINT_UNKNOWNPOINT_H_
#define SRC_MESH_ELEMENTS_POINT_UNKNOWNPOINT_H_

#define UnknownPointVTKCode 0

#include "pointelement.h"

namespace espreso {

class UnknownPoint: public PointElement
{

public:
	UnknownPoint(eslocal index): PointElement(index) {};
	Element* copy() const { return new UnknownPoint(*this); }

	eslocal vtkCode() const { return UnknownPointVTKCode; }
	size_t gaussePoints() const { ESINFO(GLOBAL_ERROR) << "Unknown point has no gausse points."; return 0; }

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



#endif /* SRC_MESH_ELEMENTS_POINT_UNKNOWNPOINT_H_ */
