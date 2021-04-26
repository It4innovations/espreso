
#include "clusterstore.h"

#include "store.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"


using namespace espreso;

ClusterStore::ClusterStore()
{
	offset = info::mpi::rank;
	size = 1;
	totalSize = info::mpi::size;
}

ClusterStore::~ClusterStore()
{

}
