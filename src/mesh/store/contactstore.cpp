
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
//: eps(0.01), groupsize(1),
//  surface(surface),
//  nodes(new NodeStore()),
//  elements(NULL),
//  enodes(NULL),
//  epointers(NULL),
//  enormals(NULL),
//  closeElements(NULL),
//  contactPairs(NULL),
//  intersections(NULL),
//  grid(NULL)
{

}

ContactStore::~ContactStore()
{
	if (localPairs) { delete localPairs; }
	if (neighPairs) { delete neighPairs; }
	if (intersections) { delete intersections; }
	if (interface) { delete interface; }
	if (planeData) { delete planeData; }
//	if (nodes != NULL) { delete nodes; }
//	if (elements != NULL) { delete elements; }
//	if (enodes != NULL) { delete enodes; }
//	if (epointers != NULL) { delete epointers; }
//	if (enormals != NULL) { delete enormals; }
//	if (closeElements != NULL) { delete closeElements; }
//	if (contactPairs != NULL) { delete contactPairs; }
//	if (intersections != NULL) { delete intersections; }
//	if (grid != NULL) { delete grid; }
//
//	for (size_t i = 0; i < gnecoords.size(); i++) {
//		delete gnecoords[i];
//	}
//	for (size_t i = 0; i < gnenodes.size(); i++) {
//		delete gnenodes[i];
//	}
//	for (size_t i = 0; i < gnepointers.size(); i++) {
//		delete gnepointers[i];
//	}
//	for (size_t i = 0; i < gneIDs.size(); i++) {
//		delete gneIDs[i];
//	}
//	for (size_t i = 0; i < gnenormals.size(); i++) {
//		delete gnenormals[i];
//	}
//	for (size_t i = 0; i < gngrid.size(); i++) {
//		delete gngrid[i];
//	}
//	for (size_t i = 0; i < gncloseElements.size(); i++) {
//		delete gncloseElements[i];
//	}
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

