
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

struct DatabaseOffset {
	esint global, local, size;

	bool operator<(const DatabaseOffset &other) { return global < other.global; }
};
typedef std::vector<DatabaseOffset> Blocks;

struct OrderedNodes {
	Blocks blocks;

	ivector<_Point<esfloat> > coordinates;
};

struct OrderedElements {
	Blocks blocks;

	ivector<Element::CODE> etype;
	ivector<esint> enodes;
};

struct OrderedFaces {
	Blocks blocks;

	esint elements;
	ivector<Element::CODE> etype;
	ivector<esint> enodes;
	ivector<esint> owner, neighbor;
};

struct OrderedValues {
	Blocks blocks;

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
	~InputMesh()
	{
		if (nodes) delete nodes;
		if (elements) delete elements;
		if (regions) delete regions;
	}
};

inline size_t size(const OrderedNodes &data)
{
	return
			data.blocks.size() * sizeof(DatabaseOffset) +
			data.coordinates.size() * sizeof(_Point<esfloat>);
}

inline size_t size(const OrderedElements &data)
{
	return
			data.blocks.size() * sizeof(DatabaseOffset) +
			data.etype.size() * sizeof(Element::CODE) +
			data.enodes.size() * sizeof(esint);
}

inline size_t size(const OrderedValues &data)
{
	return
			data.blocks.size() * sizeof(DatabaseOffset) +
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
