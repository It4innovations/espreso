
#include "fetidatastore.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

FETIDataStore::FETIDataStore()
{
	sFixPointsDistribution = {0, 0};
	iFixPointsDistribution = {0, 0};
}

FETIDataStore::~FETIDataStore()
{

}

size_t FETIDataStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(corners);

	packedSize += utils::packedSize(surfaceFixPoints);
	packedSize += utils::packedSize(sFixPointsDistribution);
	packedSize += utils::packedSize(innerFixPoints);
	packedSize += utils::packedSize(iFixPointsDistribution);

	return packedSize;
}

void FETIDataStore::packFull(char* &p) const
{
	utils::pack(corners, p);

	utils::pack(surfaceFixPoints, p);
	utils::pack(sFixPointsDistribution, p);
	utils::pack(innerFixPoints, p);
	utils::pack(iFixPointsDistribution, p);
}

void FETIDataStore::unpackFull(const char* &p)
{
	utils::unpack(corners, p);

	utils::unpack(surfaceFixPoints, p);
	utils::unpack(sFixPointsDistribution, p);
	utils::unpack(innerFixPoints, p);
	utils::unpack(iFixPointsDistribution, p);
}


