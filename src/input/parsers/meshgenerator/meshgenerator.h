
#ifndef SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_

#include "input/input.h"
#include "composition/generator.h"

#include <cstddef>

namespace espreso {

class InputGeneratorConfiguration;
class Mesh;

class MeshGenerator: public Input {
public:
	static double precision;

	MeshGenerator(InputGeneratorConfiguration &generator);
	~MeshGenerator();

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh) {}
	void finishVariables() {}
	int timeSteps() { return 0; }
	void nextTimeStep(Mesh &mesh);



protected:
	InputGeneratorConfiguration &configuration;
	Generator *generator;

	NodesDomain nodes;
	Elements elements;
	OrderedRegions regions;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_ */

