
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

/**
 * ORDERED DATA
 */

struct OrderedDataDistribution {
	esint chunk, offset, size, total;
};

struct OrderedNodesBalanced: OrderedDataDistribution {
	ivector<_Point<esfloat> > coordinates;
};

struct OrderedElementsBalanced: OrderedDataDistribution {
	ivector<Element::CODE> etype;
	ivector<esint> edist, enodes;
};

/**
 * CLUSTERD DATA
 */

struct DataDuplication {
	esint origin, duplication;
	bool operator<(const DataDuplication &other) const { if (origin == other.origin) { return duplication < other.duplication; } return origin < other.origin; }
	bool operator!=(const DataDuplication &other) const { return origin != other.origin || duplication != other.duplication; }
};

struct ClusteredDataDistribution {
	ivector<esint> offsets;
};

struct ClusteredNodes: ClusteredDataDistribution {
	ivector<_Point<esfloat> > coordinates;
};

struct MergedNodes: ClusteredNodes {
	std::vector<DataDuplication> duplication;
};

struct LinkedNodes: MergedNodes {
	std::vector<int> neighbors;
	ivector<esint> rankDistribution, rankData;
	std::unordered_map<esint, esint> g2l;
};

struct ClusteredElements: ClusteredDataDistribution {
	ivector<Element::CODE> etype;
	ivector<esint> edist, enodes;
};

struct MergedElements: ClusteredElements {
	std::vector<DataDuplication> duplication;
};

template <typename TNodes, typename TElements>
struct TemporalMesh {
	TNodes *nodes;
	TElements *elements;

	TemporalMesh()
	: nodes(new TNodes()), elements(new TElements())
	{

	}

	~TemporalMesh()
	{
		clear();
	}

	void clear()
	{
		if (nodes) { delete nodes; nodes = nullptr; }
		if (elements) { delete elements; elements = nullptr; }
	}
};

// balancing
void balance(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, int &dimension);

// clusterization
void assignBuckets(const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, const HilbertCurve<esfloat> &sfc, ivector<esint> &nbuckets, ivector<esint> &ebuckets);
void clusterize(TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, ivector<esint> &nbuckets, ivector<esint> &ebuckets, esint buckets, TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, ivector<esint> &splitters);
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, const TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, const ivector<esint> &splitters, std::vector<int> &sfcNeighbors);

// merging
void searchDuplicatedNodes(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, ClusteredNodes *clustered, MergedNodes *merged);
void searchParentAndDuplicatedElements(TemporalMesh<LinkedNodes, ClusteredElements> &linked, TemporalMesh<LinkedNodes, MergedElements> &prepared, int meshDimension);

//linking
void linkup(TemporalMesh<MergedNodes, ClusteredElements> &merged, TemporalMesh<LinkedNodes, ClusteredElements> &linked);

// filler
void fillMesh(TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
