
#include "domainstore.h"

#include "store.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"


using namespace espreso;

DomainStore::DomainStore()
: distribution({0, 0}),
  nodes(NULL),
  dual(NULL),
  localDual(NULL)
{
	size = 1;
	totalSize = info::mpi::size;
}

DomainStore::~DomainStore()
{
	if (nodes) delete nodes;
	if (dual) delete dual;
	if (localDual) delete localDual;
}

std::vector<esint> DomainStore::gatherProcDistribution()
{
	return Store::gatherDistribution(size);
}
