
#ifndef SRC_INPUT_INPUT_H_
#define SRC_INPUT_INPUT_H_

#include "input/meshbuilder.h"

#include <cstddef>
#include <string>
#include <vector>
#include <functional>

namespace espreso {

class Input {

protected:
	Input(MeshBuilder &meshData)
	: _meshData(meshData), _eregsize(1), _nregsize(1) {}

	void clip();

	void balance();
	void balanceNodes();
	void balancePermutedNodes();
	void balanceElements();
	void balancePermutedElements();

	void assignRegions(
			std::map<std::string, std::vector<esint> > &regions, std::vector<esint> &IDs,
			std::vector<esint> &distribution,
			size_t &rsize, std::vector<esint> &rbits);
	void fillRegions(std::map<std::string, std::vector<esint> > &regions, size_t &rsize, std::vector<esint> &rbits);

	void sortNodes(bool withElementNodes = false);
	void sortElements();
	void sortElements(const std::vector<esint> &permutation);

	void fillNodes();
	void fillElements();
	void fillNeighbors();

	void fillNodeRegions();
	void fillBoundaryRegions();
	void fillElementRegions();

	void reindexElementNodes();
	void reindexBoundaryNodes();

	void removeDuplicateElements();
	void searchDuplicateNodes();
	void coupleDuplicateNodes();

	void searchDuplicateNodes(std::vector<Point> &coordinates, std::vector<esint> &ids, std::function<void(esint id, esint target)> merge);

	MeshBuilder &_meshData;
	std::vector<esint> _nDistribution, _eDistribution, _etypeDistribution;

	size_t _eregsize, _nregsize;
	std::vector<esint> _eregions, _nregions;
};

}



#endif /* SRC_INPUT_INPUT_H_ */
