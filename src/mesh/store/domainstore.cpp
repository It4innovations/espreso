
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

size_t DomainStore::packedFullSize() const
{
    return
            utils::packedSize(offset) +
            utils::packedSize(size) +
            utils::packedSize(totalSize) +
            utils::packedSize(next) +
            utils::packedSize(distribution) +
            utils::packedSize(nodes) +
            utils::packedSize(elements) +
            utils::packedSize(cluster) +
            utils::packedSize(dual) +
            utils::packedSize(localDual);
}

void DomainStore::packFull(char* &p) const
{
    utils::pack(offset, p);
    utils::pack(size, p);
    utils::pack(totalSize, p);
    utils::pack(next, p);
    utils::pack(distribution, p);
    utils::pack(nodes, p);
    utils::pack(elements, p);
    utils::pack(cluster, p);
    utils::pack(dual, p);
    utils::pack(localDual, p);
}

void DomainStore::unpackFull(const char* &p)
{
    utils::unpack(offset, p);
    utils::unpack(size, p);
    utils::unpack(totalSize, p);
    utils::unpack(next, p);
    utils::unpack(distribution, p);
    utils::unpack(nodes, p);
    utils::unpack(elements, p);
    utils::unpack(cluster, p);
    utils::unpack(dual, p);
    utils::unpack(localDual, p);
}
