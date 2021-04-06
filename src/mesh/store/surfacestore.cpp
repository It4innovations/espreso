
#include "surfacestore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"
#include "esinfo/envinfo.h"

#include "mesh/mesh.h"

using namespace espreso;

SurfaceStore::SurfaceStore()
: parents(NULL),
  fID(NULL),
  contact(NULL),
  body(NULL),
  triangles(NULL),
  nodes(NULL),
  nIDs(NULL),
  coordinates(NULL),
  enodes(NULL),
  base(NULL),
  normal(NULL),
  parameters(NULL),
  nelements(NULL),
  IDs(NULL),
  neighbors(NULL),
  offset(0),
  size(0),
  totalSize(0),
  epointers(NULL)
{

}

SurfaceStore::~SurfaceStore()
{
	if (parents != NULL) { delete parents; }
	if (fID != NULL) { delete fID; }
	if (contact != NULL) { delete contact; }
	if (body != NULL) { delete body; }
	if (triangles != NULL && triangles != enodes) { delete triangles; }
	if (nodes != NULL) { delete nodes; }
	if (nIDs != NULL) { delete nIDs; }
	if (coordinates != NULL) { delete coordinates; }
	if (enodes != NULL) { delete enodes; }
	if (base != NULL) { delete base; }
	if (normal != NULL) { delete normal; }
	if (parameters != NULL) { delete parameters; }
	if (nelements != NULL) { delete nelements; }
	if (IDs != NULL) { delete IDs; }
	if (neighbors != NULL) { delete neighbors; }
	if (epointers != NULL) { delete epointers; }
}

size_t SurfaceStore::packedFullSize() const
{
	size_t packedSize = 0;

	packedSize += utils::packedSize(parents);
	packedSize += utils::packedSize(fID);
	packedSize += utils::packedSize(contact);
	packedSize += utils::packedSize(body);
	packedSize += utils::packedSize(triangles);
	packedSize += utils::packedSize(nodes);
	packedSize += utils::packedSize(nIDs);
	packedSize += utils::packedSize(coordinates);
	packedSize += utils::packedSize(enodes);
	packedSize += utils::packedSize(base);
	packedSize += utils::packedSize(normal);
	packedSize += utils::packedSize(parameters);
	packedSize += utils::packedSize(nelements);
	packedSize += utils::packedSize(IDs);
	packedSize += utils::packedSize(neighbors);

	packedSize += utils::packedSize(tdistribution);
	packedSize += utils::packedSize(edistribution);

	packedSize += utils::packedSize(offset);
	packedSize += utils::packedSize(size);
	packedSize += utils::packedSize(totalSize);
	packedSize += 1;
	if (epointers != NULL) {
		packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
	}
	packedSize += utils::packedSize(ecounters);

	return packedSize;
}

void SurfaceStore::packFull(char* &p) const
{
	utils::pack(parents, p);
	utils::pack(fID, p);
	utils::pack(contact, p);
	utils::pack(body, p);
	utils::pack(triangles, p);
	utils::pack(nodes, p);
	utils::pack(nIDs, p);
	utils::pack(coordinates, p);
	utils::pack(enodes, p);
	utils::pack(base, p);
	utils::pack(normal, p);
	utils::pack(parameters, p);
	utils::pack(nelements, p);
	utils::pack(IDs, p);
	utils::pack(neighbors, p);

	utils::pack(tdistribution, p);
	utils::pack(edistribution, p);

	utils::pack(offset, p);
	utils::pack(size, p);
	utils::pack(totalSize, p);
	utils::pack(epointers != NULL, p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());
		for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
			eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
		}
		utils::pack(eindices, p);
	}
	utils::pack(ecounters, p);
}

void SurfaceStore::unpackFull(const char* &p)
{
	utils::unpack(parents, p);
	utils::unpack(fID, p);
	utils::unpack(contact, p);
	utils::unpack(body, p);
	utils::unpack(triangles, p);
	utils::unpack(nodes, p);
	utils::unpack(nIDs, p);
	utils::unpack(coordinates, p);
	utils::unpack(enodes, p);
	utils::unpack(base, p);
	utils::unpack(normal, p);
	utils::unpack(parameters, p);
	utils::unpack(nelements, p);
	utils::unpack(IDs, p);
	utils::unpack(neighbors, p);

	utils::unpack(tdistribution, p);
	utils::unpack(edistribution, p);

	utils::unpack(offset, p);
	utils::unpack(size, p);
	utils::unpack(totalSize, p);

	bool notnull;
	utils::unpack(notnull, p);
	if (notnull) {
		std::vector<int> eindices;
		utils::unpack(eindices, p);
		if (epointers != NULL) {
			delete epointers;
		}
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, eindices.size()));
		for (size_t i = 0; i < eindices.size(); ++i) {
			epointers->datatarray()[i] = &Mesh::edata[eindices[i]];
		}
	}
	utils::unpack(ecounters, p);
}

