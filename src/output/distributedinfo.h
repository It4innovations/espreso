
#ifndef SRC_OUTPUT_DISTRIBUTEDINFO_H_
#define SRC_OUTPUT_DISTRIBUTEDINFO_H_

#include "../basis/point/point.h"
#include "meshinfo.h"

namespace espreso {

class Element;

namespace output {

struct DistributedInfo: public MeshInfo {

	DistributedInfo(const Mesh *mesh, double domainShrinkRatio, double clusterShrinkRatio, InfoMode mode = InfoMode::PREPARE);
	DistributedInfo(const Mesh *mesh, const Region* region, double domainShrinkRatio, double clusterShrinkRatio, InfoMode mode = InfoMode::PREPARE);

	MeshInfo* deriveRegion(const Region *region) const;
	MeshInfo* copyWithoutMesh() const;

	void addSettings(size_t step);
	void addSolution(const std::vector<Solution*> &solution);
	void addGeneralInfo();

	bool isShrunk() const { return _domainShrinkRatio != 1 || _clusterShrinkRatio != 1; }
	bool distributed() const { return true; }

	Point shrink(const Point &p, eslocal domain) const;

protected:
	double _domainShrinkRatio;
	double _clusterShrinkRatio;

	Point _clusterCenter;
	std::vector<Point> _domainsCenters;

	// region x domain x indices
	std::vector<std::vector<std::vector<eslocal> > > _cIndices;
	// region x domain offset
	std::vector<std::vector<size_t> > _cOffset;

private:
	void prepare(const std::vector<Element*> &region, InfoMode mode);
};

}
}



#endif /* SRC_OUTPUT_DISTRIBUTEDINFO_H_ */
