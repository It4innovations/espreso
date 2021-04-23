
#include "nodes.uniform.api.composer.h"

#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/elementstore.h"
#include "physics/system/fetisystem.h"

using namespace espreso;

NodesUniformAPIComposer::NodesUniformAPIComposer(const FETIConfiguration &configuration, int DOFs)
: NodesUniformFETIComposer(configuration, NULL, NULL, DOFs)
{

}


void NodesUniformAPIComposer::fill(FETISolverData &data)
{
	_initDOFMap();

	std::vector<esint> distribution = info::mesh->elements->gatherDomainsProcDistribution();

	data.K.initDomains(info::mesh->elements->domains.size);
	data.K.fillDecomposition(
				info::mpi::rank, info::mpi::size, info::mesh->neighbors.size(),
				distribution.data(), info::mesh->neighbors.data(), _DOFMap);
	data.f.initVectors(1);
	data.f.initDomains(DataDecomposition::DUPLICATION::SPLIT_DOMAINS, info::mesh->elements->domains.size);
	data.f.fillDecomposition(
			info::mpi::rank, info::mpi::size, info::mesh->neighbors.size(),
			distribution.data(), info::mesh->neighbors.data(), _DOFMap);
}
