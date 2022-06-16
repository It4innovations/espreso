
#ifndef SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_

#include "gridgenerator.h"

namespace espreso {

class Mesh;
struct MeshBuilder;
struct GridTowerGeneratorConfiguration;

class GridTowerGenerator: public GridGenerator {

public:
	static void generate(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh);

	static esint grids(const GridTowerGeneratorConfiguration &configuration);
	static esint gridIndex(const GridTowerGeneratorConfiguration &configuration);

protected:
	GridTowerGenerator(const GridTowerGeneratorConfiguration &configuration);
	virtual ~GridTowerGenerator() {}

	void init(const GridTowerGeneratorConfiguration &configuration);
	void nodes(MeshBuilder &mesh);
	void neighbors(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh);
	void regions(const GridTowerGeneratorConfiguration &configuration, MeshBuilder &mesh);

	esint _gridIndex;
	esint _gridNodeOffset;
};

}


#endif /* SRC_INPUT_MESHGENERATOR_COMPOSITION_GRIDTOWERGENERATOR_H_ */
