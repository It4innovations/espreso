
#ifndef SRC_INPUT_PARSERS_NGLIB_NGLIB_H_
#define SRC_INPUT_PARSERS_NGLIB_NGLIB_H_

#include "input/meshbuilder.h"

namespace espreso {

class InputConfiguration;

class NGLibGenerator: public MeshBuilder {
public:
	NGLibGenerator(const InputConfiguration &configuration);
	void load();

protected:
	const InputConfiguration &_configuration;
};

}

#endif /* SRC_INPUT_PARSERS_NGLIB_NGLIB_H_ */
