
#ifndef SRC_INPUT_PARSERS_NEPER_NEPER_H_
#define SRC_INPUT_PARSERS_NEPER_NEPER_H_

#include "input/meshbuilder.h"

namespace espreso {

class InputConfiguration;

class NeperLoader: public MeshBuilder {
public:
	NeperLoader(InputConfiguration &configuration);
	void load();

protected:
	InputConfiguration &_configuration;
};

}

#endif /* SRC_INPUT_PARSERS_NEPER_NEPER_H_ */
