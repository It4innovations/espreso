
#ifndef SRC_INPUT_BUILDERS_BUILDER_H_
#define SRC_INPUT_BUILDERS_BUILDER_H_

#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

void buildOrderedFEM(NodesBlocks &nodes, ElementsBlocks &elements, OrderedRegions &regions, Mesh &mesh);
void buildChunkedFVM(NodesBlocks &nodes, FacesBlocks &faces, OrderedRegions &regions, Mesh &mesh);
void buildOrderedFVM(NodesBlocks &nodes, FacesBlocks &faces, OrderedRegions &regions, Mesh &mesh);

void buildDecomposedFEM(NodesDomain &nodes, Elements &elements, OrderedRegions &regions, Mesh &mesh);
void buildDecomposedFVM(NodesDomain &nodes, Faces &faces, OrderedRegions &regions, Mesh &mesh);


void orderedValuesInit(VariablesBlocks &variables, Mesh &mesh);
void orderedValuesNext(VariablesBlocks &variables, Mesh &mesh);
void orderedValuesFinish(VariablesBlocks &variables);

void chunkedValuesInit(VariablesBlocks &variables, Mesh &mesh);
void chunkedValuesNext(VariablesBlocks &variables, Mesh &mesh);
void chunkedValuesFinish(VariablesBlocks &variables);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_H_ */

