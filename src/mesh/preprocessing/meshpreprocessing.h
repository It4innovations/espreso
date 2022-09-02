
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <vector>

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

// methods for 'pure' mesh
//  - elements, nodes, neighbors only
//  - no requirements on decomposition -> operations with neighbors processes can take some time
// ================================================================================================

void computeNodesDuplication(NodeStore *nodes, std::vector<int> &neighborsWithMe);

void linkNodesAndElements(ElementStore *elements, NodeStore *nodes, const std::vector<int> &neighbors);
void computeElementsFaceNeighbors(NodeStore *nodes, ElementStore *elements);
void computeElementsEdgeNeighbors(NodeStore *nodes, ElementStore *elements);

void computeElementsCenters(const NodeStore *nodes, ElementStore *elements);

int getStronglyConnectedComponents(const ElementStore *elements, std::vector<int> &component);
void computeComponentDual(const ElementStore *elements, esint coffset, esint csize, const std::vector<int> &component, const std::vector<int> &neighbors, std::vector<esint> &dualDist, std::vector<esint> &dualData);

void computeElementsClusterization(const ElementStore *elements, const NodeStore *nodes, std::vector<esint> &partition);
void exchangeElements(ElementStore* &elements, NodeStore* &nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe, const std::vector<esint> &partition);
void computeContinuousClusterization(const ElementStore *elements, const NodeStore *nodes, const std::vector<esint> &dualDist, const std::vector<esint> &dualData, esint coffset, esint csize, const std::vector<int> &component, const std::vector<int> &neighborsWithMe, std::vector<esint> &partition);

void sortNodes(NodeStore *nodes, ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions);
void computeElementDistribution(ElementStore *elements);
void computeERegionsDistribution(ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);


// methods for clustered mesh
//  - elements, nodes distribution computed
//  - operations with neighbors processes should be scalable
// ================================================================================================

void fillRegionMask(ElementStore *elements, const std::vector<ElementsRegionStore*> &elementsRegions);
void processNamelessElements(ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);

ElementStore* exchangeHalo(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors);

void computeBaseMeshInfo(ElementStore *elements, NodeStore *nodes);
void computeBodies(ElementStore *elements, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);
void fillERegionsNodesMask(const std::vector<ElementsRegionStore*> &elementsRegions, const ElementStore *elements, const std::vector<int> &neighbors, NodeStore *nodes);
void fillBRegionsNodesMask(const std::vector<BoundaryRegionStore*> &boundaryRegions, const std::vector<ContactInterfaceStore*> &contactInterfaces, const std::vector<int> &neighbors, NodeStore *nodes);
void computeERegionsNodes(const NodeStore *nodes, const ElementStore *elements, const std::vector<int> &neighbors, std::vector<ElementsRegionStore*> &elementsRegions);
void computeBRegionsNodes(const std::vector<int> &neighbors, NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeBRegionsElementsFromNodes(const NodeStore *nodes, const ElementStore *elements, const ElementStore *halo, const std::vector<ElementsRegionStore*> &elementsRegions, BoundaryRegionStore *bregion);
void computeBRegionsDistribution(NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeBRegionsParents(const NodeStore *nodes, const ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void fillOutputOffset(NodeStore *nodes, const std::vector<ElementsRegionStore*> &elementsRegions, const std::vector<BoundaryRegionStore*> &boundaryRegions, const std::vector<ContactInterfaceStore*> &contactInterfaces, const std::vector<int> &neighbors);

void computeBodiesSurface(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, SurfaceStore *surface, std::vector<int> &neighbors);
void computeWarpedNormals(SurfaceStore * surface);
void exchangeContactHalo(SurfaceStore * surface, ContactStore *contact);
void findCloseElements(ContactStore *contact);
void computeContactInterface(SurfaceStore* surface, ContactStore* contact);
void arrangeContactInterfaces(ContactStore* contact, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);

void triangularizeSurface(SurfaceStore *surface);
void triangularizeBoundary(BoundaryRegionStore *boundary);
void computeRegionsSurface(ElementStore *elements, NodeStore *nodes, ElementStore *halo, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);
void computeSurfaceElementNeighbors(NodeStore *nodes, std::vector<int> &neigbors, SurfaceStore *surface);

void computeVolumeIndices(ElementStore *elements, const NodeStore *nodes);
void computeVolumeIndicesOM(ElementStore *elements, const NodeStore *nodes);
void computeVolumeIndicesOMOpt(ElementStore *elements, const NodeStore *nodes);
void computeVolumeIndicesOMOpt2(ElementStore *elements, const NodeStore *nodes);

// methods for mesh decomposition
//  - decompose elements to domains and intervals
// ================================================================================================

void computeDecomposedDual(const ElementStore *elements, std::vector<esint> &dualDist, std::vector<esint> &dualData);
void computeElementsDecomposition(const ElementStore *elements, esint parts, std::vector<size_t> &distribution, std::vector<esint> &permutation);
void permuteElements(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces, std::vector<int> &neighbors, std::vector<size_t> &distribution, const std::vector<esint> &permutation);


// methods for clustered and decomposed mesh
//  - all data are computed
// ================================================================================================

void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints);

void computeDomainDual(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe);
void computeClustersDistribution(DomainStore *domains, ClusterStore *clusters);
void computeDomainsSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors);
void triangularizeDomainSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, SurfaceStore *domainsSurface, std::vector<int> &neighbors);

void computeNodeDomainDistribution(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<int> neighborsWithMe);
void computeLocalIndices(ElementStore *elements, DomainStore *domains);

void computeElementIntervals(const DomainStore *domains, ElementStore *elements);

void computeRegionsElementIntervals(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);
void computeRegionsBoundaryIntervals(const DomainStore *domains, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);

}
}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
