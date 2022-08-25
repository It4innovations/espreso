
#ifndef SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_GENERATOR_H_
#define SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_GENERATOR_H_

#include "input/input.h"

namespace espreso {

class Generator {
public:
	virtual ~Generator() {};

	virtual void nodes(NodesDomain &nodes) =0;
	virtual void elements(Elements &elements) =0;
	virtual void neighbors(Domain &domain) =0;
	virtual void regions() =0;
};

}



#endif /* SRC_INPUT_PARSERS_MESHGENERATOR_COMPOSITION_GENERATOR_H_ */
