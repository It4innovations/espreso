
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_

#include "triple.h"

#include <cstddef>

namespace espreso {

struct BlockGeneratorConfiguration;

struct BlockSettings {

	BlockSettings(const BlockGeneratorConfiguration &configuration);

	static size_t preferedDomains(const BlockGeneratorConfiguration &configuration);

	Triple<size_t> domains;
	Triple<size_t> elements;

	Triple<int> start, end;
};

}



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_ */
