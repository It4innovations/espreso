
#include "contactinterfacestore.h"

#include "basis/utilities/packing.h"

using namespace espreso;

ContactInterfaceStore::ContactInterfaceStore(const std::string &name, esint interfaceIndex)
: BoundaryRegionStore(name), interfaceIndex(interfaceIndex)
{

}

ContactInterfaceStore::ContactInterfaceStore(const char* &packedData)
: BoundaryRegionStore("")
{
	unpackFull(packedData);
}

size_t ContactInterfaceStore::packedFullSize() const
{
	size_t size = BoundaryRegionStore::packedFullSize();
	size += utils::packedSize(interfaceIndex);
	return size;
}

void ContactInterfaceStore::packFull(char* &p) const
{
	BoundaryRegionStore::packFull(p);
	utils::pack(interfaceIndex, p);
}

void ContactInterfaceStore::unpackFull(const char* &p)
{
	BoundaryRegionStore::unpackFull(p);
	utils::unpack(interfaceIndex, p);
}
