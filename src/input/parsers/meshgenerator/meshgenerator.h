
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
	double nextVariables(Mesh &mesh);

	int variables() { return 0; }

protected:
	InputGeneratorConfiguration &configuration;
	Generator *generator;

	NodesDomain nodes;
	Elements elements;
	OrderedRegions regions;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_ */

