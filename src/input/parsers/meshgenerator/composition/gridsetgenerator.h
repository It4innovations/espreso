
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

class Mesh;
struct MeshBuilder;
struct GridSetGeneratorConfiguration;

class GridSetGenerator: public GridGenerator {

public:
	static void generate(const GridSetGeneratorConfiguration &configuration, MeshBuilder &mesh);

	static esint grids(const GridSetGeneratorConfiguration &configuration);
	static esint gridIndex(const GridSetGeneratorConfiguration &configuration);

protected:
	GridSetGenerator(const GridSetGeneratorConfiguration &configuration);
	virtual ~GridSetGenerator() {}

	virtual void init(const GridSetGeneratorConfiguration &configuration);
	virtual void nodes(MeshBuilder &mesh);
	virtual void regions(const GridSetGeneratorConfiguration &configuration, MeshBuilder &mesh);

	esint _gridIndex;
	esint _gridNodeOffset;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDSETGENERATOR_H_ */
