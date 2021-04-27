
#ifndef SRC_INPUT_PARSERS_GMSH_GMSH_H_
#define SRC_INPUT_PARSERS_GMSH_GMSH_H_

#include "input/meshbuilder.h"

namespace espreso {

class InputConfiguration;

class GMSHGenerator: public MeshBuilder {
public:
	GMSHGenerator(const InputConfiguration &configuration);
	void load();

protected:
	const InputConfiguration &_configuration;
};

}


#endif /* SRC_INPUT_PARSERS_GMSH_GMSH_H_ */
