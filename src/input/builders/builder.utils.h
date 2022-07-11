
#ifndef SRC_INPUT_BUILDERS_BUILDER_UTILS_H_
#define SRC_INPUT_BUILDERS_BUILDER_UTILS_H_

#include "basis/sfc/hilbertcurve.h"
#include "input/input.h"
#include "mesh/mesh.h"

namespace espreso {
namespace builder {

struct ChunkDistribution {
	esint chunk, offset, size, total;
};

struct ExplicitOffset {
	ivector<esint> offsets;
};

struct DataDuplication {
	esint origin, duplicate;
	bool operator<(const DataDuplication &other) const { if (origin == other.origin) { return duplicate < other.duplicate; } return origin < other.origin; }
	bool operator!=(const DataDuplication &other) const { return origin != other.origin || duplicate != other.duplicate; }
};

struct Duplication {
	std::vector<DataDuplication> duplication;
};

struct NeighborsLinks {
	std::vector<int> neighbors;
	ivector<esint> rankDistribution, rankData;
};

/**
 * Nodes states
 *
 * 1. Ordered Nodes
 *   Elements connectivities are described by 'offsets'. Hence, nodes are stored in blocks.
 *   Nodes can be loaded arbitrary. Offsets are described by 'Blocks'.
 *
 * 2. Ordered Nodes Balanced
 *   Nodes are sorted according to 'offset' and evenly distributed across all processes in chunks.
 *
 * 3. Clustered Nodes
 *   Nodes are sorted according to 'coordinates' and evenly distributed across all processes in chunks.
 *   Original 'offsets' need to be stored explicitly.
 *
 * 4. Merged Nodes
 *   Clustered + All duplicated nodes were found.
 *
 * 5. Linked Nodes
 *   Merged + Information about neighboring 'ranks' that hold the same nodes.
 */

struct OrderedNodesBalanced: Nodes, ChunkDistribution { };
struct ClusteredNodes: Nodes, ExplicitOffset { };
struct MergedNodes: Nodes, ExplicitOffset, Duplication { };
struct LinkedNodes: Nodes, ExplicitOffset, Duplication, NeighborsLinks { };
struct NodesHolder {
	OrderedNodes ordered;
	OrderedNodesBalanced balanced;
	ClusteredNodes clustered;
	MergedNodes merged;
	LinkedNodes linked;

	NodesHolder() = default;
	NodesHolder(OrderedNodes &&ordered): ordered(std::move(ordered)) {}
};

/**
 * Elements states
 *
 * 1. Ordered Elements
 *   Elements connectivities are described by global 'offsets'.
 *   Elements can be loaded arbitrary. Offsets are described by 'Blocks'.
 *
 * 2. Ordered Elements Balanced
 *   Elements are sorted according to 'offset' and evenly distributed across all processes.
 *
 * 3. Clustered Elements
 *   Elements are sorted according to 'coordinates' and evenly distributed across all processes.
 *   Original 'offsets' need to be stored explicitly.
 *
 * 4. Merged Elements
 *   Clustered + All duplicated and parents elements were found.
 */

struct OrderedElementsBalanced: Elements, ChunkDistribution { };
struct ClusteredElements: Elements, ExplicitOffset { };
struct MergedElements: Elements, ExplicitOffset, Duplication { };
struct ElementsHolder {
	OrderedElements ordered;
	OrderedElementsBalanced balanced;
	ClusteredElements clustered;
	MergedElements merged;

	ElementsHolder() = default;
	ElementsHolder(OrderedElements &&ordered): ordered(std::move(ordered)) {}
};

struct OrderedFacesBalanced: Faces { ChunkDistribution elements, faces; };
struct FaceHolder {
	OrderedFaces ordered;
	OrderedFacesBalanced balanced;

	FaceHolder(OrderedFaces &&ordered): ordered(std::move(ordered)) {}
};

/**
 * Standard workflow for ordered database with regions.
 *
 *	  +-- NODES --+- ELEMENS -+
 *    |-----------+-----------| <-- parallel parser provides data
 * 0. | FACES     | FACES     | <-- in the case of FVM
 * 1. | ORDERED   | ORDERED   | <-- initial setting for the builder
 * 2. | BALANCED  | BALANCED  |
 * 3. |-----------------------| <-- compute SFC and assign buckets
 * 4. | CLUSTERED | CLUSTERED |
 * 5. |-----------------------| <-- approximate neighbors
 * 6. | MERGED    | CLUSTERED |
 * 7. | LINKED    | CLUSTERED |
 * 8. | LINKED    | MERGED    |
 * 9. +-----------------------+ <-- fill the mesh
 */

void swap(Nodes &n1, Nodes &n2);
void swap(Elements &e1, Elements &e2);
void swap(Faces &f1, Faces &f2);

// 0.
void trivialUpdate(OrderedFaces &ordered, OrderedFacesBalanced &balanced);
void buildElementsFromFaces(OrderedFacesBalanced &faces, OrderedElementsBalanced &elements, OrderedNodes &nodes);

// 1. -> 2.
void trivialUpdate(OrderedNodes &ordered, OrderedNodesBalanced &balanced);
void trivialUpdate(OrderedElements &ordered, OrderedElementsBalanced &balanced);
void balanceFEM(OrderedNodes &inNodes, OrderedElements &inElements, OrderedNodesBalanced &outNodes, OrderedElementsBalanced &outElements);

// 3.
void assignBuckets(OrderedNodesBalanced &nodes, OrderedElementsBalanced &elements, const HilbertCurve<esfloat> &sfc, ivector<esint> &nbuckets, ivector<esint> &ebuckets);

// 2. -> 4.
void trivialUpdate(OrderedNodesBalanced &balanced, ClusteredNodes &clustered);
void trivialUpdate(OrderedElementsBalanced &balanced, ClusteredElements &clustered);
void clusterize(OrderedNodesBalanced &inNodes, OrderedElementsBalanced &inElements, ivector<esint> &nbuckets, ivector<esint> &ebuckets, esint buckets, ClusteredNodes &outNodes, ClusteredElements &outElements, ivector<esint> &splitters);

// 4. -> 6.
void trivialUpdate(ClusteredNodes &clustered, MergedNodes &merged);
void exchangeSFCBoundaryNodes(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, const std::vector<int> &sfcNeighbors, ClusteredNodes &clustered);
void searchDuplicatedNodes(ClusteredNodes &clustered, MergedNodes &merged);

// 5.
void computeSFCNeighbors(const HilbertCurve<esfloat> &sfc, const ivector<esint> &splitters, std::vector<int> &sfcNeighbors);

// 6. -> 7.
void trivialUpdate(MergedNodes &merged, LinkedNodes &linked);
void mergeDuplicatedNodes(MergedNodes &merged); // 'linkup' has side-effect of removing duplicated nodes
void linkup(MergedNodes &merged, LinkedNodes &linked, ClusteredElements &elements);

// 7. -> 8.
void trivialUpdate(ClusteredElements &clustered, MergedElements &merged);
void mergeDuplicatedElements(ClusteredElements &clustered, MergedElements &merged, LinkedNodes &nodes, int dimension);

// 9.
void fillNodes(LinkedNodes &nodes, OrderedRegions &regions, Mesh &mesh);
void fillElements(MergedElements &elements, OrderedRegions &regions, Mesh &mesh);


// utility functions
inline size_t size(const Nodes &data)
{
	return
			data.coordinates.size() * sizeof(_Point<esfloat>);
}

inline size_t size(const Elements &data)
{
	return
			data.etype.size() * sizeof(Element::CODE) +
			data.enodes.size() * sizeof(esint);
}

inline size_t size(const Faces &data)
{
	return
			data.ftype.size() * sizeof(Element::CODE) +
			data.fnodes.size() * sizeof(esint) +
			data.owner.size() * sizeof(esint) +
			data.neighbor.size() * sizeof(esint);
}

inline size_t size(const OrderedValues &data)
{
	return
			data.blocks.size() * sizeof(DatabaseOffset) +
			data.data.size() * sizeof(esfloat);
}

inline size_t size(const OrderedRegions &data)
{
	return
			data.nodes.size() * sizeof(OrderedRegions::Region) +
			data.elements.size() * sizeof(OrderedRegions::Region);
}

}
}

#endif /* SRC_INPUT_BUILDERS_BUILDER_UTILS_H_ */
