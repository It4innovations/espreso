
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_GRID_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_GRID_H_

#include "../../loader.h"
#include "../primitives/block.h"

namespace espreso {
namespace input {

struct GridSettings {
	Triple<double> start, end;
	Triple<size_t> blocks, clusters, domains, elements;

	std::vector<bool> nonempty;
};

class Grid: public Loader {

public:
	Grid(Mesh &mesh, size_t index, size_t size);
	virtual ~Grid() { delete _block; };

	static void load(Mesh &mesh, size_t index, size_t size)
	{
		ESINFO(OVERVIEW) << "Generate grid";
		Grid grid(mesh, index, size);
		grid.fill();
	}

	virtual void points(Coordinates &coordinates);
	virtual void elements(std::vector<Element*> &elements);
	virtual void materials(std::vector<Material> &materials);
	virtual void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours);
	virtual void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

protected:
	GridSettings _generator;
	BlockGenerator* _block;
	Triple<size_t> _clusterOffset;
	Triple<size_t> _subnodes;

	size_t _index;
	size_t _size;

	std::vector<int> _cMap;
	int _clusterIndex;
};


}
}



#endif /* SRC_INPUT_GENERATOR_COMPOSITION_GRID_H_ */
