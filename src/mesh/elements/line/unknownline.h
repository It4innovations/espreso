
#ifndef SRC_MESH_ELEMENTS_LINE_UNKNOWNLINE_H_
#define SRC_MESH_ELEMENTS_LINE_UNKNOWNLINE_H_

#define UnknownLineVTKCode -1

#include "lineelement.h"

namespace espreso {

class UnknownLine: public LineElement
{
public:
	UnknownLine(std::vector<eslocal> &nodes): _nodes(nodes) {};
	Element* copy() const { return new UnknownLine(*this); }

	eslocal nCommon() const { return _nodes.size() > 4 ? 3 : 2; }
	eslocal vtkCode() const { return UnknownLineVTKCode; }

	size_t nodes() const { return _nodes.size(); }
	size_t coarseNodes() const { return _nodes.size(); }
	size_t gaussePoints() const { ESINFO(GLOBAL_ERROR) << "Unknown line has no gausse points."; return 0; }

	virtual Point faceNormal(const Element *face) const
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown line element.";
		return Point();
	}
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates)
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown line element.";
		return Point();
	}

	const std::vector<DenseMatrix>& dN() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<DenseMatrix>& N() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }
	const std::vector<double>& weighFactor() const { ESINFO(GLOBAL_ERROR) << "Unknown element has no base functions"; exit(EXIT_FAILURE); }

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
	eslocal* indices() { return _nodes.data(); }
	const eslocal* indices() const { return _nodes.data(); }

private:
	std::vector<eslocal> &_nodes;
};


}

#endif /* SRC_MESH_ELEMENTS_LINE_UNKNOWNLINE_H_ */
