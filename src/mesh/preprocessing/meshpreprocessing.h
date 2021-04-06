
#ifndef SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_
#define SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_

#include <cstddef>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct Element;
struct RegionStore;
struct BoundaryRegionStore;
struct SurfaceStore;
template <typename TEBoundaries, typename TEData> class serializededata;

namespace mesh {

void linkNodesAndElements();
void exchangeHalo();
void exchangeElements(const std::vector<esint> &partition);

void computeElementsFaceNeighbors();
void computeElementsEdgeNeighbors();
void computeElementsCenters();
void computeDecomposedDual(std::vector<esint> &dualDist, std::vector<esint> &dualData);
void computeDomainDual();

void reclusterize();
void partitiate(esint parts, bool uniformDecomposition);
void computeBodies();

void arrangeNodes();
void arrangeElements();
void arrangeElementsRegions();
void arrangeBoundaryRegions();

void computeNodeDomainDistribution();
void computeLocalIndices();
void computeDomainsSurface();
void triangularizeDomainSurface();

void triangularizeSurface(SurfaceStore *surface);
void triangularizeBoundary(BoundaryRegionStore *boundary);

void computeRegionsSurface();
void computeBodiesSurface();
void computeSurfaceLocations();
void computeSurfaceElementNeighbors(SurfaceStore *surface);
void computeContactNormals();
void fireNormals();
void findCloseElements();
//void clip(std::vector<Point> &p, std::vector<Point> &q, std::vector<std::vector<Point> > &res);
//void dummyClip(std::vector<Point> &p, std::vector<Point> &q, std::vector<std::vector<Point> > &res);

void computeBoundaryElementsFromNodes(BoundaryRegionStore *bregion, int elementDimension);

void linkNodesAndElements(
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		std::vector<size_t> &edistribution,
		bool sortedIDs);

void computeElementsNeighbors(
		serializededata<esint, esint>* &nelements,
		serializededata<esint, esint>* &eneighbors,
		serializededata<esint, esint> *enodes,
		serializededata<esint, esint> *eIDs,
		serializededata<esint, Element*> *epointers,
		std::vector<size_t> &edistribution,
		std::function<serializededata<int, int>*(Element*)> across,
		bool insertNeighSize,
		bool sortedIDs);

esint getSFCDecomposition(std::vector<esint> &partition);
esint callParallelDecomposer(std::vector<esint> &eframes, std::vector<esint> &eneighbors, std::vector<esint> &partition);
void permuteElements(const std::vector<esint> &permutation, const std::vector<size_t> &distribution);
void arrangeElementsPermutation(std::vector<esint> &permutation);
void computeBoundaryNodes(std::vector<esint> &externalBoundary, std::vector<esint> &internalBoundary);
void fillRegionMask();

void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints);

void synchronizeRegionNodes(std::vector<RegionStore*> &regions);
void computeNodeInfo(std::vector<RegionStore*> &regions);

}
}



#endif /* SRC_MESH_PREPROCESSING_MESHPREPROCESSING_H_ */
