
#ifndef SRC_OUTPUT_COLLECTEDINFO_H_
#define SRC_OUTPUT_COLLECTEDINFO_H_

#include "../basis/point/point.h"
#include "meshinfo.h"

namespace espreso {

class Element;

namespace output {

struct CollectedInfo: public MeshInfo {

	CollectedInfo(const Mesh *mesh, InfoMode mode = InfoMode::PREPARE);
	CollectedInfo(const Mesh *mesh, const Region* region, InfoMode mode = InfoMode::PREPARE);

	MeshInfo* deriveRegion(const Region *region) const;
	MeshInfo* copyWithoutMesh() const;

	void addSettings(size_t step);
	void addSolution(const std::vector<Solution*> &solution);
	void addGeneralInfo();

	bool isShrunk() const { return false; }
	bool distributed() const { return false; }

	Point shrink(const Point &p, eslocal domain) const { return p; }

protected:
	std::vector<esglobal> _globalIDs; // ID0end, ID2end, ID3end, ...
	std::vector<esglobal> _globalIDsMap; // ID0position0, ID0position1, ..., ID1position0, ...
	std::vector<esglobal> _globalIDsMultiplicity;

	// region x indices
	std::vector<std::vector<eslocal> > _cIndices;

private:
	void prepare(const std::vector<Element*> &region);

};

}
}


#endif /* SRC_OUTPUT_COLLECTEDINFO_H_ */
