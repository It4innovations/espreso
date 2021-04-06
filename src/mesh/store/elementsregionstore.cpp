
#include "elementsregionstore.h"
#include "surfacestore.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

using namespace espreso;

ElementsRegionStore::ElementsRegionStore(const std::string &name)
: RegionStore(name),
  elements(NULL),
  uniqueElements(NULL),
  surface(new SurfaceStore())
{

}

ElementsRegionStore::ElementsRegionStore(const char* &packedData)
: ElementsRegionStore("")
{
	unpackFull(packedData);
}

ElementsRegionStore::~ElementsRegionStore()
{
	if (uniqueElements != NULL && uniqueElements != elements) {
		delete uniqueElements;
	}
	if (elements != NULL) { delete elements; }
	delete surface;
}

size_t ElementsRegionStore::packedFullSize() const
{
	size_t packedSize = RegionStore::packedFullSize();

	packedSize += utils::packedSize(elements);
	packedSize += utils::packedSize(uniqueElements);

	packedSize += utils::packedSize(eintervals);
	packedSize += utils::packedSize(ueintervals);

	return packedSize;
}

void ElementsRegionStore::packFull(char* &p) const
{
	RegionStore::packFull(p);
	utils::pack(elements, p);
	utils::pack(uniqueElements, p);
	utils::pack(eintervals, p);
	utils::pack(ueintervals, p);
}

void ElementsRegionStore::unpackFull(const char* &p)
{
	RegionStore::unpackFull(p);
	utils::unpack(elements, p);
	utils::unpack(uniqueElements, p);
	utils::unpack(eintervals, p);
	utils::unpack(ueintervals, p);
}

size_t ElementsRegionStore::packedSize() const
{
	if (elements == NULL) {
		return 0;
	}
	return
			RegionStore::packedSize() +
			elements->packedSize() +
			utils::packedSize(eintervals);
}

void ElementsRegionStore::pack(char* &p) const
{
	if (elements == NULL) {
		return;
	}
	RegionStore::pack(p);
	elements->pack(p);
	utils::pack(eintervals, p);
}

void ElementsRegionStore::unpack(const char* &p)
{
	RegionStore::unpack(p);
	if (elements == NULL) {
		elements = new serializededata<esint, esint>(1, tarray<esint>(1, 0));
	}
	elements->unpack(p);
	utils::unpack(eintervals, p);
}

