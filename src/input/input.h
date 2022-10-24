
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

	virtual void initVariables(Mesh &mesh) = 0;
	virtual void finishVariables() = 0;
	virtual int variables() = 0;
	virtual void nextVariables(Mesh &mesh) = 0;

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
	esint elements;
};

struct Variable {
	int dimension;
	std::string name;
	std::vector<esfloat> data;

	Variable(int dimenstion, const std::string &name): dimension(dimenstion), name(name) {}
};

struct VariableLoader;

struct Variables {
	std::vector<Variable> nodes, elements;

	VariableLoader *loader; // filled by builder
};

// Databases with ordered nodes and elements
struct DatabaseOffset {
	esint global, local, size;
};

struct IregularVector {
	ivector<esint> distribution, data;
};

struct Blocks {
	std::vector<DatabaseOffset> blocks;
};

struct Domain {
	ivector<esint> ids;
	IregularVector ranks;
	std::vector<int> neighbors;
};

struct NodesDomain: Nodes, Domain { };
struct NodesBlocks: Nodes, Blocks { };
struct ElementsBlocks: Elements, Blocks { };
struct FacesBlocks: Faces { Blocks edist, fdist; };
struct VariablesBlocks: Variables { Blocks edist, ndist; };

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
