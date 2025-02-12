
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

size_t ClusterStore::packedFullSize() const
{
    return
            utils::packedSize(offset) +
            utils::packedSize(size) +
            utils::packedSize(totalSize) +
            utils::packedSize(next);
}

void ClusterStore::packFull(char* &p) const
{
    utils::pack(offset, p);
    utils::pack(size, p);
    utils::pack(totalSize, p);
    utils::pack(next, p);
}

void ClusterStore::unpackFull(const char* &p)
{
    utils::unpack(offset, p);
    utils::unpack(size, p);
    utils::unpack(totalSize, p);
    utils::unpack(next, p);
}
