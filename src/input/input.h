
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "basis/containers/allocators.h"
#include "basis/containers/point.h"
#include "basis/io/inputfile.h"
#include "mesh/element.h"

#include <string>
#include <vector>

namespace espreso {

struct InputConfiguration;
class Mesh;

class Input {
public:
	virtual void load(const InputConfiguration &configuration) = 0;
	virtual void build(Mesh &mesh) = 0;
	virtual void variables(Mesh &mesh) = 0;

	virtual ~Input() {}
};

// offset -> index to vector of data
struct DatabaseOffset {
	esint global;
	esint local, size;

	bool operator<(const DatabaseOffset &other) { return global < other.global; }
};

struct OrderedUniqueNodes {
	ivector<_Point<esfloat> > coordinates;
};

struct OrderedUniqueElements {
	ivector<Element::CODE> etype;
	ivector<esint> enodes;
};

struct OrderedUniqueFaces {
	esint elements;
	ivector<Element::CODE> etype;
	ivector<esint> enodes;
	ivector<esint> owner, neighbor;
};

struct OrderedNodes {
	std::vector<DatabaseOffset> offsets;
	ivector<_Point<esfloat> > coordinates;
};

struct OrderedElements {
	std::vector<DatabaseOffset> offsets;
	ivector<Element::CODE> etype;
	ivector<esint> enodes;
};

struct OrderedValues {
	std::vector<DatabaseOffset> offsets;
	ivector<esfloat> data;
};

struct OrderedRegions {
	struct Region {
		std::string name;
		int dimension;
		esint start, end;
	};

	std::vector<Region> nodes, elements;
};

template <typename TNodes, typename TElements, typename TRegions>
struct InputMesh {
	TNodes *nodes;
	TElements *elements;
	TRegions *regions;

	InputMesh(): nodes(new TNodes()), elements(new TElements()), regions(new TRegions()) { }
};

inline size_t size(const OrderedUniqueNodes &data)
{
	return data.coordinates.size() * sizeof(_Point<esfloat>);
}

inline size_t size(const OrderedUniqueElements &data)
{
	return data.etype.size() * sizeof(Element::CODE) + data.enodes.size() * sizeof(esint);
}

inline size_t size(const OrderedUniqueFaces &data)
{
	return
			data.etype.size() * sizeof(Element::CODE) +
			data.enodes.size() * sizeof(esint) +
			data.owner.size() * sizeof(esint) +
			data.neighbor.size() * sizeof(esint);
}

inline size_t size(const OrderedNodes &data)
{
	return
			data.offsets.size() * sizeof(DatabaseOffset) +
			data.coordinates.size() * sizeof(_Point<esfloat>);
}

inline size_t size(const OrderedElements &data)
{
	return
			data.offsets.size() * sizeof(DatabaseOffset) +
			data.etype.size() * sizeof(Element::CODE) +
			data.enodes.size() * sizeof(esint);
}

inline size_t size(const OrderedValues &data)
{
	return
			data.offsets.size() * sizeof(DatabaseOffset) +
			data.data.size() * sizeof(esfloat);
}

inline size_t size(const OrderedRegions &data)
{
	return
			data.nodes.size() * sizeof(OrderedRegions::Region) +
			data.elements.size() * sizeof(OrderedRegions::Region);
}

}

#endif /* SRC_INPUT_INPUT_H_ */
