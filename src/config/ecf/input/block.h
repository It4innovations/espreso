
#ifndef SRC_CONFIG_ECF_INPUT_BLOCK_H_
#define SRC_CONFIG_ECF_INPUT_BLOCK_H_

#include "generatorelements.h"
#include "config/description.h"
#include "config/holders/expression.h"

#include <string>

namespace espreso {

struct BlockGeneratorConfiguration: public ECFDescription {

	GENERATOR_ELEMENT_TYPE element_type;

	ECFExpression start_x, start_y, start_z;
	ECFExpression length_x, length_y, length_z;

	std::string projection_x, projection_y, projection_z;
	std::string rotation_x, rotation_y, rotation_z;

	size_t domains_x, domains_y, domains_z;
	size_t elements_x, elements_y, elements_z;

	BlockGeneratorConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_BLOCK_H_ */
