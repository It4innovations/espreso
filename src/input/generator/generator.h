
#ifndef SRC_INPUT_GENERATOR_FACTORY_H_
#define SRC_INPUT_GENERATOR_FACTORY_H_

#include <cstddef>

namespace espreso {

struct ESPRESOGenerator;
class Mesh;

namespace input {

class Generator {

	static void generate(const ESPRESOGenerator &configuration, Mesh &mesh, size_t index, size_t size);
};
}
}



#endif /* SRC_INPUT_GENERATOR_FACTORY_H_ */
