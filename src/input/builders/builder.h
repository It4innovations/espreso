
#ifndef SRC_INPUT_BUILDERS_BUILDER_H_
#define SRC_INPUT_BUILDERS_BUILDER_H_

#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

void build(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_H_ */

