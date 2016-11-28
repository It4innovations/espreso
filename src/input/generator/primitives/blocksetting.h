
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_

#include "esbasis.h"

namespace espreso {
namespace input {

struct BlockSetting {
	Triple<size_t> domains;
	Triple<size_t> elements;

	Triple<double> start, end;
	Triple<Expression> projection = Triple<Expression>(Expression("x", { "x", "y", "z" }), Expression("y", { "x", "y", "z" }), Expression("z", { "x", "y", "z" }));
	Triple<Expression> rotation = Triple<Expression>(Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }), Expression("0", { "x", "y", "z" }));
};

}
}



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_ */
