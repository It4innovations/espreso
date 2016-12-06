
#ifndef SRC_INPUT_GENERATOR_GENERATOR_H_
#define SRC_INPUT_GENERATOR_GENERATOR_H_

#include <cstddef>

namespace espreso {

struct GlobalConfiguration;
class Mesh;

namespace input {

struct Generator {

	static void generate(const GlobalConfiguration &configuration, Mesh &mesh, size_t index, size_t size);
};
}
}



#endif /* SRC_INPUT_GENERATOR_GENERATOR_H_ */
