
#ifndef SRC_INPUT_BUILDERS_BUILDER_UTILS_H_
#define SRC_INPUT_BUILDERS_BUILDER_UTILS_H_

#include "basis/sfc/hilbertcurve.h"
#include "input/input.h"
#include "mesh/mesh.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace espreso {
namespace builder {

struct PackedData {
	std::vector<esint> distribution, data;
};

struct DataDuplication {
	esint origin, duplication;

	bool operator<(const DataDuplication &other) const { if (origin == other.origin) { return duplication < other.duplication; } return origin < other.origin; }
};

struct OrderedMesh: public OrderedMeshDatabase {
	int dimension;
	esint nchunk, noffset, nsize, ntotal;
	esint echunk, eoffset, esize, etotal;

	std::vector<esint> edist;
	PackedData ndata, edata;
};

struct ClusteredMesh: public OrderedMeshDatabase {
	int dimension;
	size_t buckets;
	std::vector<esint, initless_allocator<esint> > nbuckets, ebuckets;
	std::vector<esint> splitters;
	std::vector<esint> noffsets, eoffsets;

	std::vector<esint> edist;
	std::vector<int> neighbors;

	std::unordered_map<esint, esint> g2l;
	std::vector<DataDuplication> nduplication, eduplication;
	std::vector<esint> rankDistribution, rankData;
};

// balancing
void balance(OrderedMeshDatabase &database, OrderedMesh &mesh);

// clusterization
void assignBuckets(const OrderedMesh &mesh, const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);
void clusterize(OrderedMesh &mesh, ClusteredMesh &clustered);
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);

// merging
void searchDuplicatedNodes(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);
void searchParentAndDuplicatedElements(ClusteredMesh &linked);

//linking
void linkup(ClusteredMesh &clustered, ClusteredMesh &linked);

// filler
void fillNodes(ClusteredMesh &linked, Mesh &mesh);
void fillElements(ClusteredMesh &linked, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
