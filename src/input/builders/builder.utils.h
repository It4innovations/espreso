
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

struct OrderedFacesBalanced: OrderedElementsBalanced {
	ivector<esint> owner, neighbor, foffset;
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

inline size_t size(const OrderedNodesBalanced &data)
{
	return data.coordinates.size() * sizeof(_Point<esfloat>);
}

inline size_t size(const OrderedElementsBalanced &data)
{
	return
			data.etype.size() * sizeof(Element::CODE) +
			data.edist.size() * sizeof(esint) +
			data.enodes.size() * sizeof(esint);
}

template <typename TNodes, typename TElements>
struct TemporalMesh {
	TNodes *nodes;
	TElements *elements;

	TemporalMesh()
	: nodes(new TNodes()), elements(new TElements()),
	  _nodes(nodes), _elements(elements)
	{

	}

	template <typename TNodesOther, typename TElementsOther>
	TemporalMesh(const TemporalMesh<TNodesOther, TElementsOther> &other)
	: nodes(dynamic_cast<TNodes*>(other.nodes)),
	  elements(dynamic_cast<TElements*>(other.elements)),
	  _nodes(nullptr), _elements(nullptr)
	{

	}

	~TemporalMesh()
	{
		if (_nodes) { delete _nodes; }
		if (_elements) { delete _elements; }
	}

private:
	// instance holder that should be deleted at the end
	TNodes *_nodes;
	TElements *_elements;
};

template <typename TNodes, typename TElements>
struct TemporalSequentialMesh: TemporalMesh<TNodes, TElements> {};

// sequential
void initializeSequentialFEM(const InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, const TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, int &dimension);
void initializeSequentialFVM(const InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, const TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> &clustered);

// balancing
void balanceFEM(const InputMesh<OrderedNodes, OrderedElements, OrderedRegions> &input, const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, int &dimension);
void balanceFVM(const InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> &input, const TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> &ordered);

// clusterization
void assignBuckets(const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, const HilbertCurve<esfloat> &sfc, ivector<esint> &nbuckets, ivector<esint> &ebuckets);
void clusterize(const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered, ivector<esint> &nbuckets, ivector<esint> &ebuckets, esint buckets, const TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, ivector<esint> &splitters);
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, std::vector<int> &sfcNeighbors);

// merging
void searchDuplicatedNodes(const TemporalSequentialMesh<ClusteredNodes, ClusteredElements> &clustered, const TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged);
void searchDuplicatedElements(const TemporalSequentialMesh<MergedNodes, ClusteredElements> &merged, const TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, int meshDimension);

void searchDuplicatedNodesWithSFC(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, const TemporalMesh<ClusteredNodes, ClusteredElements> &clustered, const TemporalMesh<MergedNodes, ClusteredElements> &merged);
void searchParentAndDuplicatedElements(const TemporalMesh<LinkedNodes, ClusteredElements> &linked, const TemporalMesh<LinkedNodes, MergedElements> &prepared, int meshDimension);
void reindexToLocal(const TemporalMesh<LinkedNodes, MergedElements> &linked);

void buildElementsFromFaces(const TemporalSequentialMesh<MergedNodes, OrderedFacesBalanced> &clustered, const TemporalSequentialMesh<MergedNodes, MergedElements> &prepared);
void buildElementsFromFaces(const TemporalMesh<OrderedNodesBalanced, OrderedFacesBalanced> &grouped, const TemporalMesh<OrderedNodesBalanced, OrderedElementsBalanced> &ordered);

//linking
void linkup(const TemporalMesh<MergedNodes, ClusteredElements> &merged, const TemporalMesh<LinkedNodes, ClusteredElements> &linked);

// filler
void fillMesh(const TemporalMesh<LinkedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh);
void fillSequentialMesh(const TemporalSequentialMesh<MergedNodes, MergedElements> &prepared, OrderedRegions &regions, Mesh &mesh);

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
