
#ifndef SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_
#define SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_

#include "input/meshbuilder.h"

#include <cstddef>

namespace espreso {

class InputGeneratorConfiguration;
class Mesh;

class MeshGenerator: public MeshBuilder {
public:
	static double precision;

	MeshGenerator(const InputGeneratorConfiguration &configuration);
	void load();

protected:
	const InputGeneratorConfiguration &_configuration;
};

}



#endif /* SRC_INPUT_MESHGENERATOR_MESHGENERATOR_H_ */
