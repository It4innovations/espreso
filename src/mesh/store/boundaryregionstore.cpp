
#include "boundaryregionstore.h"

#include "esinfo/envinfo.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/packing.h"

#include "mesh/mesh.h"

using namespace espreso;

BoundaryRegionStore::BoundaryRegionStore(const std::string &name)
: RegionStore(name),
  originalDimension(0),
  dimension(0),
  area(0),
  procNodes(NULL),
  triangles(NULL),
  epointers(NULL),
  emembership(NULL)
{

}

BoundaryRegionStore::BoundaryRegionStore(const char* &packedData)
: BoundaryRegionStore("")
{
	unpackFull(packedData);
}

BoundaryRegionStore::~BoundaryRegionStore()
{
	if (procNodes != NULL) { delete procNodes; }
	if (triangles != NULL && triangles != procNodes) { delete triangles; }

	if (epointers != NULL) { delete epointers; }
	if (emembership != NULL) { delete emembership; }
}

void BoundaryRegionStore::permute(const std::vector<esint> &permutation, const std::vector<size_t> &distribution)
{
	this->distribution = distribution;

	if (procNodes != NULL) {
		procNodes->permute(permutation, distribution);
	}

	if (epointers != NULL) {
		epointers->permute(permutation, distribution);
	}
	if (emembership != NULL) {
		emembership->permute(permutation, distribution);
	}

	eintervals.clear();
}

size_t BoundaryRegionStore::packedFullSize() const
{
	size_t packedSize = RegionStore::packedFullSize();

	packedSize += utils::packedSize(originalDimension);
	packedSize += utils::packedSize(dimension);
	packedSize += utils::packedSize(area);

	packedSize += utils::packedSize(distribution);
	packedSize += utils::packedSize(procNodes);
	packedSize += utils::packedSize(triangles);

	packedSize += 1;
	if (epointers != NULL) {
		packedSize += sizeof(size_t) + epointers->datatarray().size() * sizeof(int);
	}
	packedSize += 1;
	if (emembership != NULL) {
		packedSize += sizeof(size_t) + emembership->datatarray().size() * sizeof(esint);
	}

	packedSize += utils::packedSize(eintervals);
	packedSize += utils::packedSize(eintervalsDistribution);

	return packedSize;
}

void BoundaryRegionStore::packFull(char* &p) const
{
	RegionStore::packFull(p);
	utils::pack(originalDimension, p);
	utils::pack(dimension, p);
	utils::pack(area, p);

	utils::pack(distribution, p);
	utils::pack(procNodes, p);
	utils::pack(triangles, p);

	utils::pack(epointers != NULL, p);
	if (epointers != NULL) {
		std::vector<int> eindices;
		eindices.reserve(epointers->datatarray().size());
		for (size_t i = 0; i < epointers->datatarray().size(); ++i) {
			eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
		}
		utils::pack(eindices, p);
	}
	utils::pack(emembership, p);

	utils::pack(eintervals, p);
	utils::pack(eintervalsDistribution, p);
}

void BoundaryRegionStore::unpackFull(const char* &p)
{
	RegionStore::unpackFull(p);
	utils::unpack(originalDimension, p);
	utils::unpack(dimension, p);
	utils::unpack(area, p);

	utils::unpack(distribution, p);
	utils::unpack(procNodes, p);
	utils::unpack(triangles, p);

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
	utils::unpack(emembership, p);

	utils::unpack(eintervals, p);
	utils::unpack(eintervalsDistribution, p);
}

size_t BoundaryRegionStore::packedSize() const
{
	return
			RegionStore::packedSize() +
			utils::packedSize(originalDimension) +
			utils::packedSize(dimension) +
			procNodes->packedSize() +
			sizeof(size_t) + epointers->datatarray().size() * sizeof(int) +
			utils::packedSize(emembership) +
			utils::packedSize(eintervals);
}

void BoundaryRegionStore::pack(char* &p) const
{
	RegionStore::pack(p);
	utils::pack(originalDimension, p);
	utils::pack(dimension, p);
	procNodes->pack(p);
	std::vector<int> eindices;
	eindices.reserve(epointers->datatarray().size());

	size_t threads = info::env::OMP_NUM_THREADS;
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = epointers->datatarray().distribution()[t]; i < epointers->datatarray().distribution()[t + 1]; ++i) {
			eindices.push_back(static_cast<int>(epointers->datatarray()[i]->code));
		}
	}
	utils::pack(eindices, p);
	utils::pack(emembership, p);
	utils::pack(eintervals, p);
}

void BoundaryRegionStore::unpack(const char* &p)
{
	RegionStore::unpack(p);
	utils::unpack(originalDimension, p);
	utils::unpack(dimension, p);
	if (procNodes == NULL) {
		procNodes = new serializededata<esint, esint>(tarray<esint>(1, 0), tarray<esint>(1, 0));
		epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, 0));
	}
	procNodes->unpack(p);
	std::vector<int> eindices;
	utils::unpack(eindices, p);
	if (epointers != NULL) {
		delete epointers;
	}
	epointers = new serializededata<esint, Element*>(1, tarray<Element*>(1, eindices.size()));
	for (size_t i = 0; i < eindices.size(); ++i) {
		epointers->datatarray()[i] = &Mesh::edata[eindices[i]];
	}
	utils::unpack(emembership, p);
	utils::unpack(eintervals, p);
}
