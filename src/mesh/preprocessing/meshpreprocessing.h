
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct Element;
struct ElementStore;
struct NodeStore;
struct RegionStore;
struct DomainStore;
struct BoundaryRegionStore;
struct ElementsRegionStore;
struct ContactInterfaceStore;
struct SurfaceStore;
struct ContactStore;
struct BodyStore;
struct ClusterStore;
struct ElementsInterval;
struct ElementsDistributionInfo;
template <typename TEBoundaries, typename TEData> class serializededata;

namespace mesh {

void triangularizeSurface(SurfaceStore *surface);
void triangularizeBoundary(BoundaryRegionStore *boundary);

void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints);


/// new methods

inline int bitMastSize(size_t elements)
{
	size_t size = 8 * sizeof(esint);
	return elements / size + (elements % size ? 1 : 0);
}

ElementStore* exchangeHalo(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors);

void computeNodesDuplication(NodeStore *nodes, std::vector<int> &neighborsWithMe);

void linkNodesAndElements(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors);
void linkNodesAndElements(
		NodeStore *nodes,
		std::vector<int> &neighbors,
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		std::vector<size_t> &edistribution,
		bool sortedIDs);

void computeElementsNeighbors(
		NodeStore *nodes,
		std::vector<int> &neighbors,
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint>* &eneighbors,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		serializededata<esint, Element*> *epointers,
		std::vector<size_t> &edistribution,
		std::function<serializededata<int, int>*(Element*)> across,
		bool insertNeighSize,
		bool sortedIDs);

void computeElementsCenters(ElementStore *elements, NodeStore *nodes);

void computeElementsFaceNeighbors(NodeStore *nodes, ElementStore *elements, std::vector<int> &neighbors);
void computeElementsEdgeNeighbors(NodeStore *nodes, ElementStore *elements, std::vector<int> &neighbors);

void computeDecomposedDual(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors, std::vector<esint> &dualDist, std::vector<esint> &dualData);
void computeDomainDual(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe);
void computeDomainsSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors);
void triangularizeDomainSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors);

void computeNodeDomainDistribution(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<int> neighborsWithMe);
void computeLocalIndices(ElementStore *elements, DomainStore *domains);

void computeRegionsSurface(ElementStore *elements, NodeStore *nodes, ElementStore *halo, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);
void computeSurfaceElementNeighbors(NodeStore *nodes, std::vector<int> &neigbors, SurfaceStore *surface);

void computeBoundaryNodes(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<esint> &externalBoundary, std::vector<esint> &internalBoundary);

void sortNodes(NodeStore *nodes, ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions);

void computeElementIntervals(const DomainStore *domains, ElementStore *elements);
void computeElementDistribution(const std::vector<ElementsInterval> &eintervals, ElementsDistributionInfo &distribution);

void computeRegionsElementNodes(const NodeStore *nodes, const ElementStore *elements, const std::vector<int> &neighbors, std::vector<ElementsRegionStore*> &elementsRegions);
void computeRegionsElementIntervals(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);
void computeRegionsElementDistribution(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);

void computeRegionsBoundaryNodes(const std::vector<int> &neighbors, NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeRegionsBoundaryParents(const NodeStore *nodes, const ElementStore *elements, const DomainStore *domains, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeRegionsBoundaryDistribution(std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeRegionsBoundaryElementsFromNodes(const NodeStore *nodes, const ElementStore *elements, const ElementStore *halo, const std::vector<ElementsRegionStore*> &elementsRegions, BoundaryRegionStore *bregion);

void computeBodies(NodeStore *nodes, ElementStore *elements, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);

void computeBodiesSurface(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, SurfaceStore *surface, std::vector<int> &neighbors);
void computeWarpedNormals(SurfaceStore * surface);
void exchangeContactHalo(SurfaceStore * surface, ContactStore *contact);
void findCloseElements(ContactStore *contact);
void computeContactInterface(SurfaceStore* surface, ContactStore* contact);
void arrangeContactInterfaces(ContactStore* contact, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);

void synchronizeRegionNodes(const NodeStore *nodes, const std::vector<int> &neighbors, std::vector<RegionStore*> &regions);
void computeNodeInfo(const NodeStore *nodes, const std::vector<int> &neighbors, std::vector<RegionStore*> &regions);

void fillRegionMask(const ElementsDistributionInfo &distribution, const std::vector<ElementsRegionStore*> &elementsRegions, serializededata<esint, esint>* &mask);

esint getSFCDecomposition(ElementStore *elements, NodeStore *nodes, std::vector<esint> &partition);
esint callParallelDecomposer(ElementStore *elements, NodeStore *nodes, std::vector<esint> &eframes, std::vector<esint> &eneighbors, std::vector<esint> &partition);
void reclusterize(ElementStore *elements, NodeStore *nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe);
void partitiate(ElementStore *elements, NodeStore *nodes, ClusterStore *clusters, DomainStore *domains, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, esint parts, bool uniformDecomposition);

void exchangeElements(ElementStore *elements, NodeStore *nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe, const std::vector<esint> &partition);
void permuteElements(ElementStore *elements, NodeStore *nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors, const std::vector<esint> &permutation, const std::vector<size_t> &distribution);

}
}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
