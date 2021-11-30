
#ifndef SRC_INPUT_BUILDERS_BUILDER_UTILS_H_
#define SRC_INPUT_BUILDERS_BUILDER_UTILS_H_

#include "basis/sfc/hilbertcurve.h"
#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

struct PackedData {
	std::vector<esint> distribution, data;
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
	std::vector<size_t> typeDistribution;
	std::vector<int> neighbors;
};

// balancing
void balance(OrderedMeshDatabase &database, OrderedMesh &mesh);

// clusterization
void assignBuckets(const OrderedMesh &mesh, const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);
void clusterize(OrderedMesh &mesh, ClusteredMesh &clustered);
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);

// merging
void mergeDuplicatedNodes(const HilbertCurve<esfloat> &sfc, ClusteredMesh &clustered);

// sorting
void groupElementTypes(ClusteredMesh &clustered);

//linking
void linkup(ClusteredMesh &clustered);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
