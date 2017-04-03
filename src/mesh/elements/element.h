#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <algorithm>

#include "../../basis/logging/logging.h"

#define __WITHOUT_DOF__ -1
#define __HAS_DOF__ 1

namespace espreso {

class Point;
class Region;
class Coordinates;
class Mesh;
class DenseMatrix;
enum class Property;

class Element
{
	friend class Mesh;

public:
	enum class Type {
		POINT = 0,
		LINE = 1,
		PLANE = 2,
		VOLUME = 3,
	};

	enum class ElementPointType {
		GAUSSE_POINT,
		VERTEX_POINT
	};

	enum Params {
		MATERIAL,
		CONSTANT,
		COORDINATES,
		BODY,
		PARAMS_SIZE
	};

	inline static bool match(const eslocal *indices, eslocal x, eslocal y)
	{
		return indices[x] == indices[y];
	}

	void store(std::ofstream& os, const Coordinates &coordinates, size_t part);
	friend std::ostream& operator<<(std::ostream& os, const Element &e);

	bool operator<(const Element& other) const
	{
		if (coarseNodes() != other.coarseNodes()) {
			return coarseNodes() < other.coarseNodes();
		}

		std::vector<eslocal> e1(indices(), indices() + coarseNodes());
		std::vector<eslocal> e2(other.indices(), other.indices() + other.coarseNodes());
		std::sort(e1.begin(), e1.end());
		std::sort(e2.begin(), e2.end());
		return e1 < e2;
	}

	bool operator==(const Element& other) const
	{
		if (coarseNodes() != other.coarseNodes()) {
			return false;
		}
		return std::is_permutation(indices(), indices() + coarseNodes(), other.indices());
	}

	virtual ~Element() {};

	virtual const std::vector<DenseMatrix>& facedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;
	virtual const std::vector<DenseMatrix>& faceN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;
	virtual const std::vector<DenseMatrix>& edgedN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;
	virtual const std::vector<DenseMatrix>& edgeN(size_t index, ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;

	virtual const std::vector<DenseMatrix>& dN(ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;
	virtual const std::vector<DenseMatrix>& N(ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;
	virtual const std::vector<double>& weighFactor(ElementPointType type = ElementPointType::GAUSSE_POINT) const = 0;

	virtual const std::vector<Property>& elementDOFs() const = 0;
	virtual const std::vector<Property>& faceDOFs() const = 0;
	virtual const std::vector<Property>& edgeDOFs() const = 0;
	virtual const std::vector<Property>& pointDOFs() const = 0;
	virtual const std::vector<Property>& midPointDOFs() const = 0;

	virtual eslocal nCommon() const = 0;
	virtual eslocal vtkCode() const = 0;

	virtual size_t filledFaces() const = 0;
	virtual size_t filledEdges() const = 0;
	virtual size_t faces() const = 0;
	virtual size_t edges() const = 0;
	virtual size_t nodes() const = 0;
	virtual size_t coarseNodes() const = 0;
	virtual size_t gaussePoints() const = 0;

	virtual const std::vector<eslocal>& faceNodes(size_t index) const =0;
	virtual const std::vector<eslocal>& edgeNodes(size_t index) const =0;

	virtual Element* face(size_t index) const = 0;
	virtual Element* edge(size_t index) const = 0;
	eslocal& node(size_t index) { return indices()[index]; }
	const eslocal& node(size_t index) const { return indices()[index]; }
	virtual eslocal* indices() = 0;
	virtual const eslocal* indices() const = 0;

	virtual eslocal param(Params param) const =0;
	virtual void setParam(Params param, eslocal value) =0;
	virtual size_t params() const =0;

	virtual Type type() const =0;

	std::vector<Region*>& regions() { return _regions; }
	const std::vector<Region*>& regions() const { return _regions; }

	std::vector<Element*>& parentElements() { return _parentElements; }
	const std::vector<Element*>& parentElements() const { return _parentElements; }

	std::vector<Element*>& parentFaces() { return _parentFaces; }
	const std::vector<Element*>& parentFaces() const { return _parentFaces; }

	std::vector<Element*>& parentEdges() { return _parentEdges; }
	const std::vector<Element*>& parentEdges() const { return _parentEdges; }

	std::vector<eslocal>& domains() { return _domains; }
	const std::vector<eslocal>& domains() const { return _domains; }

	std::vector<eslocal>& clusters() { return _clusters; }
	const std::vector<eslocal>& clusters() const { return _clusters; }

	virtual std::vector<eslocal>& DOFsIndices() { return _DOFsIndices; }
	virtual const std::vector<eslocal>& DOFsIndices() const { return _DOFsIndices; }

	std::vector<eslocal>& DOFsDomainsCounters() { return _DOFsDomainsCounters; }
	const std::vector<eslocal>& DOFsDomainsCounters() const { return _DOFsDomainsCounters; }

	std::vector<eslocal>& clusterOffsets() { return _clusterOffsets; }

	virtual const std::vector<double>& stiffnessMatrix() const { ESINFO(GLOBAL_ERROR) << "Stiffness matrix of an element is not set."; exit(0); }

	size_t DOFOffset(eslocal domain, eslocal DOFIndex) const
	{
		for (size_t i = 0; i < _DOFsIndices.size() / _domains.size(); i++) {
			if (DOFIndex == this->DOFIndex(domain, i)) {
				return i;
			}
		}
		ESINFO(ERROR) << "Element has not requested DOF index";
		return -1;
	}

	eslocal DOFIndex(eslocal domain, size_t DOFoffset) const
	{
		auto it = std::lower_bound(_domains.begin(), _domains.end(), domain);
		auto DOFs = _DOFsIndices.size() / _domains.size();
		return _DOFsIndices[DOFs * (it - _domains.begin()) + DOFoffset];
	}

	eslocal DOFCounter(eslocal cluster, size_t DOFoffset) const
	{
		auto it = std::lower_bound(_clusters.begin(), _clusters.end(), cluster);
		auto DOFs = _DOFsDomainsCounters.size() / _clusters.size();
		return _DOFsDomainsCounters[DOFs * (it - _clusters.begin()) + DOFoffset];
	}

	eslocal clusterOffset(eslocal cluster) const
	{
		return _clusterOffsets[lower_bound(_clusters.begin(), _clusters.end(), cluster)  -_clusters.begin()];
	}

	size_t numberOfGlobalDomainsWithDOF(size_t DOFoffset) const
	{
		size_t n = 0;
		for (size_t c = 0; c < _clusters.size(); c++) {
			n += _DOFsDomainsCounters[c * _DOFsDomainsCounters.size() / _clusters.size() + DOFoffset];
		}
		return n;
	}

	size_t numberOfLocalDomainsWithDOF(size_t DOFoffset) const
	{
		size_t n = 0;
		for (size_t d = 0; d < _domains.size(); d++) {
			if (_DOFsIndices[d * _DOFsIndices.size() / _domains.size() + DOFoffset] != __WITHOUT_DOF__) {
				n++;
			}
		}
		return n;
	}

	void numberOfGlobalDomains(size_t size) { _domainsCounters = size; }
	size_t numberOfGlobalDomains() const { return _domainsCounters; }

	bool inDomain(eslocal domain)
	{
		auto it = std::lower_bound(_domains.begin(), _domains.end(), domain);
		return it != _domains.end() && *it == domain;
	}

	void addParent(Element* parent) { _parentElements.push_back(parent); }
	virtual void addFace(Element* face) = 0;
	virtual void addEdge(Element* edge) = 0;

	virtual Element* addFace(const std::vector<eslocal> &nodes) = 0;

	void rotateOutside(const Element* parent, const Coordinates &coordinates, Point &normal) const;

	bool hasProperty(Property property, size_t step) const;
	double sumProperty(Property property, eslocal node, size_t step, double defaultValue) const;
	double getProperty(Property property, eslocal node, size_t step, double defaultValue) const;

protected:
	virtual Element* copy() const =0;
	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;

	virtual void setFace(size_t index, Element* face) = 0;
	virtual void setEdge(size_t index, Element* edge) = 0;

	virtual size_t fillFaces() = 0;
	virtual size_t fillEdges() = 0;


	template <class TEdge>
	Element* addUniqueEdge(const eslocal *line, size_t filled, size_t coarseSize)
	{
		for (size_t i = 0; i < filled; i++) {
			if (std::is_permutation(edge(i)->indices(), edge(i)->indices() + coarseSize, line)) {
				return this->edge(i);
			}
		}
		Element *e = new TEdge(line);
		addEdge(e);
		e->addParent(this);
		return e;
	}

	template <class TFace>
	Element* addUniqueFace(const eslocal *face, size_t filled, size_t coarseSize)
	{
		for (size_t i = 0; i < filled; i++) {
			if (std::is_permutation(this->face(i)->indices(), this->face(i)->indices() + coarseSize, face)) {
				return this->face(i);
			}
		}
		Element *f = new TFace(face);
		addFace(f);
		f->addParent(this);
		return f;
	}

	std::vector<Region*> _regions;
	std::vector<Element*> _parentElements;
	std::vector<Element*> _parentFaces;
	std::vector<Element*> _parentEdges;
	std::vector<eslocal> _domains;
	std::vector<eslocal> _clusters;
	std::vector<eslocal> _DOFsIndices;
	std::vector<eslocal> _DOFsDomainsCounters;
	std::vector<eslocal> _clusterOffsets;
	size_t _domainsCounters;
};

inline std::ostream& operator<<(std::ostream& os, const Element &e)
{
	for (size_t i = 0; i < e.nodes(); i++) {
		os << e.node(i) << " ";
	}
	return os;
}

}


#endif /* ELEMENT_H_ */
