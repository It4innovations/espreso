
#include "decomposition.direct.h"

#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void DecompositionDirect::update(std::vector<int> &neighbors)
{
    this->neighbors = neighbors;
    neighDOF.resize(neighbors.size() + 1, begin); // the last is my offset
    std::vector<esint> dBuffer = { begin };
    if (!Communication::gatherUniformNeighbors(dBuffer, neighDOF, neighbors)) {
        eslog::internalFailure("cannot exchange matrix decomposition info.\n");
    }
}


