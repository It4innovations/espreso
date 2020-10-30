
#include "contactstore.h"
#include "nodestore.h"

#include "basis/containers/serializededata.h"

using namespace espreso;

ContactStore::ContactStore(SurfaceStore *surface)
: surface(surface),
  localPairs(NULL),
  neighPairs(NULL),
  intersections(NULL),
  interface(NULL),
  planeData(NULL)
{

}

ContactStore::~ContactStore()
{
	if (localPairs) { delete localPairs; }
	if (neighPairs) { delete neighPairs; }
	if (intersections) { delete intersections; }
	if (interface) { delete interface; }
	if (planeData) { delete planeData; }
}

size_t ContactStore::packedFullSize() const
{
	size_t packedSize = 0;
	return packedSize;
}

void ContactStore::packFull(char* &p) const
{

}

void ContactStore::unpackFull(const char* &p)
{

}

