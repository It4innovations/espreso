
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_

#include "triple.h"
#include "basis/expression/expression.h"

namespace espreso {

struct BlockGeneratorConfiguration;

struct BlockSettings {

	BlockSettings(const BlockGeneratorConfiguration &configuration);

	static size_t preferedDomains(const BlockGeneratorConfiguration &configuration);

	Triple<size_t> domains;
	Triple<size_t> elements;

	Triple<int> start, end;
	Triple<Expression> projection = Triple<Expression>(Expression("x", { "x", "y", "z" }), Expression("y", { "x", "y", "z" }), Expression("z", { "x", "y", "z" }));
	Triple<Expression> rotation = Triple<Expression>(Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }));
};

}



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_ */
