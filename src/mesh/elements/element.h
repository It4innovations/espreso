#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "esbasis.h"
#include "../settings/setting.h"

namespace espreso {

class Mesh;

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

	void store(std::ofstream& os, const Coordinates &coordinates, size_t part)
	{
		eslocal value = vtkCode(), pSize = params(), p;
		os.write(reinterpret_cast<const char *>(&value), sizeof(eslocal));
		for (size_t n = 0; n < nodes(); n++) {
			eslocal index = coordinates.localIndex(node(n), part);
			os.write(reinterpret_cast<const char *>(&index), sizeof(eslocal));
		}
		os.write(reinterpret_cast<const char *>(&pSize), sizeof(eslocal));
		if (pSize) {
			p = param(Element::MATERIAL);
			os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
			p = param(Element::CONSTANT);
			os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
			p = param(Element::COORDINATES);
			os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
			p = param(Element::BODY);
			os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		}
	}
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

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>& N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	virtual const std::vector<Property>& elementDOFs() const = 0;
	virtual const std::vector<Property>& faceDOFs() const = 0;
	virtual const std::vector<Property>& edgeDOFs() const = 0;
	virtual const std::vector<Property>& pointDOFs() const = 0;
	virtual const std::vector<Property>& midPointDOFs() const = 0;

	virtual eslocal nCommon() const = 0;
	virtual eslocal vtkCode() const = 0;

	virtual size_t faces() const = 0;
	virtual size_t edges() const = 0;
	virtual size_t nodes() const = 0;
	virtual size_t coarseNodes() const = 0;
	virtual size_t gaussePoints() const = 0;

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

	Settings& settings() { return _settings; }
	const Settings& settings() const { return _settings; }

	void addSettings(Property property, Evaluator* evaluator) { return _settings[property].push_back(evaluator); }
	const std::vector<Evaluator*>& settings(Property property) const { return _settings[property]; }

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

	eslocal DOFIndex(eslocal domain, size_t DOFIndex) const
	{
		auto it = std::lower_bound(_domains.begin(), _domains.end(), domain);
		auto DOFs = _DOFsIndices.size() / _domains.size();
		return _DOFsIndices[DOFs * (it - _domains.begin()) + DOFIndex];
	}

	eslocal DOFCounter(eslocal cluster, size_t DOFIndex) const
	{
		auto it = std::lower_bound(_clusters.begin(), _clusters.end(), cluster);
		auto DOFs = _DOFsDomainsCounters.size() / _clusters.size();
		return _DOFsDomainsCounters[DOFs * (it - _clusters.begin()) + DOFIndex];
	}

	eslocal clusterOffset(eslocal cluster) const
	{
		return _clusterOffsets[lower_bound(_clusters.begin(), _clusters.end(), cluster)  -_clusters.begin()];
	}

	size_t numberOfGlobalDomainsWithDOF(size_t index) const
	{
		size_t n = 0;
		for (size_t c = 0; c < _clusters.size(); c++) {
			n += _DOFsDomainsCounters[c * _DOFsDomainsCounters.size() / _clusters.size() + index];
		}
		return n;
	}

	size_t numberOfLocalDomainsWithDOF(size_t index) const
	{
		size_t n = 0;
		for (size_t d = 0; d < _domains.size(); d++) {
			if (_DOFsIndices[d * _DOFsIndices.size() / _domains.size() + index] != -1) {
				n++;
			}
		}
		return n;
	}

	bool inDomain(eslocal domain)
	{
		auto it = std::lower_bound(_domains.begin(), _domains.end(), domain);
		return it != _domains.end() && *it == domain;
	}

	virtual void addFace(Element* face) = 0;
	virtual void addEdge(Element* edge) = 0;

	void rotateOutside(const Element* parent, const Coordinates &coordinates, Point &normal) const
	{
		Point eMid(0, 0, 0), mid(0, 0, 0);
		for (size_t i = 0; i < parent->coarseNodes(); i++) {
			eMid += coordinates[parent->node(i)];
		}
		eMid /= parent->coarseNodes();

		for (size_t i = 0; i < coarseNodes(); i++) {
			mid += coordinates[node(i)];
		}
		mid /= coarseNodes();

		Point outside = mid - eMid;
		if (outside.x * normal.x + outside.y * normal.y + outside.z * normal.z > 0) {
			normal.flip();
		}
	}

protected:
	virtual Element* copy() const =0;
	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;

	virtual void setFace(size_t index, Element* face) = 0;
	virtual void setEdge(size_t index, Element* edge) = 0;

	virtual size_t fillFaces() = 0;
	virtual size_t fillEdges() = 0;


	template <class TEdge>
	void addUniqueEdge(eslocal *line, size_t filled)
	{
		for (size_t i = 0; i < filled; i++) {
			if (std::is_permutation(edge(i)->indices(), edge(i)->indices() + 2, line)) {
				return;
			}
		}
		addEdge(new TEdge(line));
	}

	template <class TFace>
	void addUniqueFace(eslocal *face, size_t filled, size_t coarseSize)
	{
		for (size_t i = 0; i < filled; i++) {
			if (std::is_permutation(this->face(i)->indices(), this->face(i)->indices() + coarseSize, face)) {
				return;
			}
		}
		addFace(new TFace(face));
	}

	Settings _settings;
	std::vector<Element*> _parentElements;
	std::vector<Element*> _parentFaces;
	std::vector<Element*> _parentEdges;
	std::vector<eslocal> _domains;
	std::vector<eslocal> _clusters;
	std::vector<eslocal> _DOFsIndices;
	std::vector<eslocal> _DOFsDomainsCounters;
	std::vector<eslocal> _clusterOffsets;
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
