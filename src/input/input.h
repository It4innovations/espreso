
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

struct Nodes {
	ivector<_Point<esfloat> > coordinates;
};

struct Elements {
	ivector<Element::CODE> etype;
	ivector<esint> enodes;
};

struct Faces {
	ivector<Element::CODE> ftype;
	ivector<esint> fnodes;
	ivector<esint> owner, neighbor;
};

struct Values {
	ivector<esfloat> data;
};

// Databases with ordered nodes and elements
struct DatabaseOffset {
	esint global, local, size;
};

struct Blocks {
	std::vector<DatabaseOffset> blocks;
};

struct OrderedNodes: Nodes, Blocks { };
struct OrderedElements: Elements, Blocks { };
struct OrderedFaces: Faces { Blocks elements, faces; };
struct OrderedValues: Values, Blocks { };

struct OrderedRegions {
	struct Region {
		std::string name;
		int dimension;
		esint start, end;
	};

	std::vector<Region> nodes, elements;
};

}

#endif /* SRC_INPUT_INPUT_H_ */
