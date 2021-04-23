
#include "domainstore.h"
#include "store.h"

#include "esinfo/mpiinfo.h"

using namespace espreso;

DomainStore::DomainStore()
: distribution({0, 0})
{
	size = 1;
	totalSize = info::mpi::size;
}

DomainStore::~DomainStore()
{

}

std::vector<esint> DomainStore::gatherProcDistribution()
{
	return Store::gatherDistribution(size);
}
