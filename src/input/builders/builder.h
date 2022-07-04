
#ifndef SRC_INPUT_BUILDERS_BUILDER_H_
#define SRC_INPUT_BUILDERS_BUILDER_H_

#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

//void buildOrderedFEM(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh);
//void buildOrderedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh);
//void buildDecomposedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh);

void buildOrderedFEM(OrderedNodes &nodes, OrderedElements &elements, OrderedRegions &regions, Mesh &mesh);
void buildDecomposedFVM(OrderedNodes &nodes, OrderedFaces &elements, OrderedRegions &regions, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_H_ */

