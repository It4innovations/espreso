
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_

namespace espreso {
namespace input {

struct BlockSetting {
	Triple<size_t> domains;
	Triple<size_t> elements;

	Triple<double> start, end;
};

}
}



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_BLOCKSETTING_H_ */
