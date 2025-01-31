
#include "decomposition.feti.h"

#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

void DecompositionFETI::update(std::vector<int> &neighbors)
{
    neighDomain.resize(neighbors.size() + 1);
    auto ddist = info::mesh->domains->gatherProcDistribution(); // remove this
    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        neighDomain[n] = ddist[info::mesh->neighbors[n]];
    }
    neighDomain.back() = dbegin;
    DecompositionDirect::update(neighbors);
}
