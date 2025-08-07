
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct MaterialConfiguration;
struct Element;
struct ElementStore;
struct NodeStore;
struct RegionStore;
struct DomainStore;
struct DomainSurfaceStore;
struct BoundaryRegionStore;
struct ElementsRegionStore;
struct ContactInterfaceStore;
struct SurfaceStore;
struct ContactStore;
struct BodyStore;
struct ClusterStore;
struct ElementsInterval;
struct ElementsDistributionInfo;
struct FETIDataStore;
template <typename TEBoundaries, typename TEData> class serializededata;

namespace mesh {

// methods for 'pure' mesh
//  - elements, nodes, neighbors only
//  - no requirements on decomposition -> operations with neighbors processes can take some time
// ================================================================================================

void computeNodesDuplication(NodeStore *nodes, std::vector<int> &neighborsWithMe);

void linkNodesAndElements(ElementStore *elements, NodeStore *nodes, const std::vector<int> &neighbors);
void linkNodesAndElements(
        const NodeStore *nodes,
        const std::vector<int> &neighbors,
        serializededata<esint, esint>* &nelements,
        serializededata<esint, esint> *enodes,
        serializededata<esint, esint> *eIDs,
        bool sortedIDs);

void computeElementsFaceNeighbors(NodeStore *nodes, ElementStore *elements, const std::vector<int> &neighbors);
void computeElementsEdgeNeighbors(NodeStore *nodes, ElementStore *elements, const std::vector<int> &neighbors);
void computeElementsNeighbors(
        NodeStore *nodes,
        const std::vector<int> &neighbors,
        serializededata<esint, esint>* &nelements,
        serializededata<esint, esint>* &eneighbors,
        serializededata<esint, esint> *enodes,
        serializededata<esint, esint> *eIDs,
        serializededata<esint, Element*> *epointers,
        std::function<serializededata<int, int>*(Element*)> across,
        bool insertNeighSize,
        bool sortedIDs);

void computeElementsCenters(const NodeStore *nodes, ElementStore *elements);

int getStronglyConnectedComponents(const ElementStore *elements, std::vector<int> &component);
void computeComponentDual(const ElementStore *elements, esint coffset, esint csize, const std::vector<int> &component, const std::vector<int> &neighbors, std::vector<esint> &dualDist, std::vector<esint> &dualData);

void computeElementsClusterization(const ElementStore *elements, BodyStore *bodies, const NodeStore *nodes, std::vector<esint> &partition, std::vector<int> &neighbors);
void exchangeElements(ElementStore* &elements, NodeStore* &nodes, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe, const std::vector<esint> &partition);
void computeContinuousClusterization(const ElementStore *elements, const NodeStore *nodes, const std::vector<esint> &dualDist, const std::vector<esint> &dualData, esint coffset, esint csize, const std::vector<int> &component, const std::vector<int> &neighborsWithMe, std::vector<esint> &partition);

void sortNodes(NodeStore *nodes, ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions);
void computeElementDistribution(ElementStore *elements);
void computeRegionsElementDistribution(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);


// methods for clustered mesh
//  - elements, nodes distribution computed
//  - operations with neighbors processes should be scalable
// ================================================================================================

void fillRegionMask(ElementStore *elements, const std::vector<ElementsRegionStore*> &elementsRegions);
void processNamelessElements(ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);

ElementStore* exchangeHalo(ElementStore *elements, NodeStore *nodes, std::vector<int> &neighbors);

void computeBodies(ElementStore *elements, BodyStore *bodies, std::vector<int> &neighbors);
void computeBodies(ElementStore *elements, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);
void computeRegionsElementNodes(const NodeStore *nodes, const ElementStore *elements, const std::vector<int> &neighbors, std::vector<ElementsRegionStore*> &elementsRegions);
void computeRegionsBoundaryNodes(const std::vector<int> &neighbors, NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeRegionsBoundaryElementsFromNodes(const NodeStore *nodes, const ElementStore *elements, const ElementStore *halo, const std::vector<ElementsRegionStore*> &elementsRegions, BoundaryRegionStore *bregion);
void computeRegionsBoundaryDistribution(NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);
void computeRegionsBoundaryParents(const NodeStore *nodes, const ElementStore *elements, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);

void computeBodiesSurface(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, SurfaceStore *surface, std::vector<int> &neighbors);
void computeSurfaceNodeNormals(NodeStore *nodes, SurfaceStore * surface, const std::vector<int> &neighbors, const double* displacement = nullptr);
void computeWarpedNormals(SurfaceStore * surface, const double* displacement = nullptr);
void exchangeContactHalo(SurfaceStore * surface, ContactStore *contact, const double* displacement = nullptr);
void findCloseElements(ContactStore *contact, const double* displacement = nullptr);
void computeContactInterface(SurfaceStore* surface, ContactStore* contact, const double *displacement = nullptr);
void arrangeContactInterfaces(NodeStore *nodes, ContactStore* contact, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<ContactInterfaceStore*> &contactInterfaces, const double* displacement = nullptr);

void triangularizeSurface(SurfaceStore *surface);
void triangularizeBoundary(BoundaryRegionStore *boundary);
void computeRegionsSurface(ElementStore *elements, NodeStore *nodes, ElementStore *halo, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<int> &neighbors);
void computeSurfaceElementNeighbors(NodeStore *nodes, std::vector<int> &neigbors, SurfaceStore *surface);


// methods for mesh decomposition
//  - decompose elements to domains and intervals
// ================================================================================================

void computeDecomposedDual(const ElementStore *elements, std::vector<esint> &dualDist, std::vector<esint> &dualData);
void computeElementsDecomposition(const ElementStore *elements, esint parts, std::vector<size_t> &distribution, std::vector<esint> &permutation);
void permuteElements(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces, std::vector<int> &neighbors, std::vector<size_t> &distribution, const std::vector<esint> &permutation);


// methods for clustered and decomposed mesh
//  - all data are computed
// ================================================================================================

void computeDomainDual(NodeStore *nodes, ElementStore *elements, DomainStore *domains, std::vector<int> &neighbors, std::vector<int> &neighborsWithMe);
void computeClustersDistribution(DomainStore *domains, ClusterStore *clusters);
void computeDomainsSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, DomainSurfaceStore *domainsSurface, std::vector<int> &neighbors);
//void triangularizeDomainSurface(NodeStore *nodes, ElementStore *elements, DomainStore *domains, DomainSurfaceStore *domainsSurface, std::vector<int> &neighbors);

void computeNodeDomainDistribution(ElementStore *elements, NodeStore *nodes, DomainStore *domains, std::vector<int> neighborsWithMe);
void computeLocalIndices(ElementStore *elements, DomainStore *domains);

void computeElementIntervals(const DomainStore *domains, ElementStore *elements);
void setMaterialsToRegions(ElementStore *elements, const std::vector<ElementsRegionStore*> &elementsRegions, const std::vector<MaterialConfiguration*> &materials, const std::map<std::string, std::string> &material_set);

void computeRegionsElementIntervals(const ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions);
void computeRegionsBoundaryIntervals(const ElementStore *elements, const DomainStore *domains, std::vector<BoundaryRegionStore*> &boundaryRegions, std::vector<ContactInterfaceStore*> &contactInterfaces);

}
}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
