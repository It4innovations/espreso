
#include "bodystore.h"

#include "store.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"


using namespace espreso;

BodyStore::BodyStore()
{
	size = 1;
	totalSize = info::mpi::size;
}

BodyStore::~BodyStore()
{

}

size_t BodyStore::packedFullSize() const
{
	return
			utils::packedSize(offset) +
			utils::packedSize(size) +
			utils::packedSize(totalSize) +
			utils::packedSize(next);
}

void BodyStore::packFull(char* &p) const
{
	utils::pack(offset, p);
	utils::pack(size, p);
	utils::pack(totalSize, p);
	utils::pack(next, p);
}

void BodyStore::unpackFull(const char* &p)
{
	utils::unpack(offset, p);
	utils::unpack(size, p);
	utils::unpack(totalSize, p);
	utils::unpack(next, p);
}
