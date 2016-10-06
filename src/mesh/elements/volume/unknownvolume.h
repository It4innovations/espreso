
#ifndef SRC_MESH_ELEMENTS_VOLUME_UNKNOWNVOLUME_H_
#define SRC_MESH_ELEMENTS_VOLUME_UNKNOWNVOLUME_H_

#define UnknownVolumeVTKCode -3

#include "volumeelement.h"

namespace espreso {

class UnknownVolume: public VolumeElement
{
public:
	UnknownVolume(std::vector<eslocal> &nodes): _nodes(nodes) {};
	Element* copy() const { return new UnknownVolume(*this); }

	eslocal nCommon() const { return _nodes.size() > 8 ? 4 : 3; }
	eslocal vtkCode() const { return UnknownVolumeVTKCode; }
	eslocal param(Params param) const { ESINFO(GLOBAL_ERROR) << "Call param of unknown volume element."; return -1; }
	void setParam(Params param, eslocal value) { ESINFO(GLOBAL_ERROR) << "Set param of unknown volume element."; }
	size_t params() const { return 0; }

	size_t faces() const { return _faces.size(); }
	size_t edges() const { return _edges.size(); }
	size_t nodes() const { return _nodes.size(); }
	size_t coarseNodes() const { return _nodes.size(); }
	size_t gaussePoints() const { ESINFO(GLOBAL_ERROR) << "Unknown volume has no gausse points."; return 0; }

	virtual Point faceNormal(const Element *face) const
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown volume element.";
		return Point();
	}
	virtual Point edgeNormal(const Element *edge, const Coordinates &coordinates) const
	{
		ESINFO(GLOBAL_ERROR) << "Call normal of unknown volume element.";
		return Point();
	}

	virtual Element* face(size_t index) const { return _faces[index]; }
	virtual Element* edge(size_t index) const { return _edges[index]; }

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
		ESINFO(GLOBAL_ERROR) << "Call neighbour of unknown volume element node.";
		return std::vector<eslocal>();
	}
	eslocal* indices() { return _nodes.data(); }
	const eslocal* indices() const { return _nodes.data(); }

	void setFace(size_t index, Element* face) { _faces[index] = face; }
	void setEdge(size_t index, Element* edge) { _edges[index] = edge; }
	void setFace(Element* face) { _faces.push_back(face); }
	void setEdge(Element* edge) { _edges.push_back(edge); }

	void fillFaces();
	void fillEdges() { ESINFO(GLOBAL_ERROR) << "Unknown volume element cannot fill edges."; }

private:
	std::vector<eslocal> &_nodes;
};


}


#endif /* SRC_MESH_ELEMENTS_VOLUME_UNKNOWNVOLUME_H_ */
