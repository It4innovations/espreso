
#ifndef SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_
#define SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_

#define UnknownPlaneVTKCode -2

#include "planeelement.h"

namespace espreso {

class UnknownPlane: public PlaneElement
{
public:
	UnknownPlane(const std::vector<Element*> &nodes, std::vector<eslocal> &indices, std::vector<eslocal> &DOFsIndices, std::vector<double> &stiffnessMatrix)
	: _nodes(nodes), _indices(indices), _stiffnessMatrix(stiffnessMatrix) {  _DOFsIndices = DOFsIndices; }
	UnknownPlane(const std::vector<Element*> &nodes, std::vector<eslocal> &indices, std::vector<double> &stiffnessMatrix)
	: _nodes(nodes), _indices(indices), _stiffnessMatrix(stiffnessMatrix) {};
	Element* copy() const { return new UnknownPlane(*this); }

	eslocal nCommon() const { return _indices.size() > 4 ? 3 : 2; }
	eslocal vtkCode() const { return UnknownPlaneVTKCode; }

	size_t edges() const { return _edges.size(); }
	size_t nodes() const { return _indices.size(); }
	size_t coarseNodes() const { return _indices.size(); }
	size_t gaussePoints() const { ESINFO(GLOBAL_ERROR) << "Unknown plane has no gausse points."; return 0; }

	const std::vector<double>& stiffnessMatrix() const { return _stiffnessMatrix; }

	const std::vector<eslocal>& edgeNodes(size_t index) const
	{
		// It is impossible to compute it
		static std::vector<eslocal> _edgeNodes;
		return _edgeNodes;
	}

	const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Unknown plane has no base functions for edge."; exit(1); }
	const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(ERROR) << "Unknown plane has no base functions for edge."; exit(1); }

	const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }

	const std::vector<Property>& elementDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& faceDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& edgeDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& pointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }
	const std::vector<Property>& midPointDOFs() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no DOFs"; exit(EXIT_FAILURE); }

protected:
	std::vector<eslocal> getNeighbours(size_t nodeIndex) const
	{
		ESINFO(GLOBAL_ERROR) << "Call neighbour of unknown plane element node.";
		return std::vector<eslocal>();
	}
	eslocal* indices() { return _indices.data(); }
	const eslocal* indices() const { return _indices.data(); }

	size_t fillEdges();

private:
	const std::vector<Element*> &_nodes;
	std::vector<eslocal> &_indices;
	std::vector<double> &_stiffnessMatrix;
	std::vector<std::vector<eslocal> > _edgeNodes;
};


}

#endif /* SRC_MESH_ELEMENTS_PLANE_UNKNOWNPLANE_H_ */
