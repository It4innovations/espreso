
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
	esint origin, duplicate;
	bool operator<(const DataDuplication &other) const { if (origin == other.origin) { return duplicate < other.duplicate; } return origin < other.origin; }
	bool operator!=(const DataDuplication &other) const { return origin != other.origin || duplicate != other.duplicate; }
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

template <typename TNodes, typename TElements>
struct TemporalSequentialMesh: TemporalMesh<TNodes, TElements> {};

// sequential
void initialize(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, int &dimension);
void initialize(InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedUniqueFacesRegions> &input, TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, int &dimension);


// balancing
void balance(InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, int &dimension);

// clusterization
void assignBuckets(const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, const HilbertCurve<esfloat> &sfc, ivector<esint> &nbuckets, ivector<esint> &ebuckets);
void clusterize(TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, ivector<esint> &nbuckets, ivector<esint> &ebuckets, esint buckets, TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, ivector<esint> &splitters);
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, std::vector<int> &sfcNeighbors);

// merging
void searchDuplicatedNodes(TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged);
void searchDuplicatedElements(TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged, TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, int meshDimension);

void searchDuplicatedNodes(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, TemporalMesh<MergedNodes, ClusteredElements> &merged);
void searchParentAndDuplicatedElements(TemporalMesh<LinkedNodes, ClusteredElements> &linked, TemporalMesh<LinkedNodes, MergedElements> &prepared, int meshDimension);

void buildElementsFromFaces();

//linking
void linkup(TemporalMesh<MergedNodes, ClusteredElements> &merged, TemporalMesh<LinkedNodes, ClusteredElements> &linked);

// filler
void fillMesh(TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh);
void fillMesh(TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
