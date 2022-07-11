
#ifndef SRC_INPUT_BUILDERS_BUILDER_H_
#define SRC_INPUT_BUILDERS_BUILDER_H_

#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

//void buildOrderedFEM(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, Mesh &mesh);
//void buildOrderedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh);
//void buildDecomposedFVM(InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> &input, Mesh &mesh);

void buildOrderedFEM(NodesBlocks &nodes, ElementsBlocks &elements, OrderedRegions &regions, Mesh &mesh);
void buildChunkedFVM(NodesBlocks &nodes, FacesBlocks &faces, OrderedRegions &regions, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_H_ */

