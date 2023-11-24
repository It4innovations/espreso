
#ifndef SRC_INPUT_PARSERS_NETGEN_NETGEN_H_
#define SRC_INPUT_PARSERS_NETGEN_NETGEN_H_

#include "input/meshbuilder.h"

namespace espreso {

struct InputConfiguration;

class NetgenNeutralLoader: public MeshBuilder {
public:
	NetgenNeutralLoader(const InputConfiguration &configuration);
	void load();

protected:
	const InputConfiguration &_configuration;
};

}

#endif /* SRC_INPUT_PARSERS_NETGEN_NETGEN_H_ */
