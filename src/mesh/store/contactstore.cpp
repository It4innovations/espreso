
#include "contactstore.h"
#include "nodestore.h"

#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"

using namespace espreso;

ContactStore::ContactStore()
: pairs(NULL),
  intersections(NULL),
  sparseSide(NULL),
  denseSide(NULL),
  planeCoordinates(NULL)
{

}

ContactStore::~ContactStore()
{
    // the last surface is the local surface
    for (size_t i = 0; i + 1 < surfaces.size(); ++i) {
        delete surfaces[i];
    }
    if (pairs) { delete pairs; }
    if (intersections) { delete intersections; }
    if (sparseSide) { delete sparseSide; }
    if (denseSide) { delete denseSide; }
    if (planeCoordinates) { delete planeCoordinates; }
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

