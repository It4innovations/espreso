
#ifndef SRC_OUTPUT_DISTRIBUTEDINFO_H_
#define SRC_OUTPUT_DISTRIBUTEDINFO_H_

#include "regioninfo.h"
#include "../basis/point/point.h"

namespace espreso {
namespace output {

struct DistributedInfo: public RegionInfo {

	DistributedInfo(const Mesh *mesh, size_t body, double domainShrinkRatio, double clusterShrinkRatio);
	DistributedInfo(const Mesh *mesh, const Region* region, double domainShrinkRatio, double clusterShrinkRatio);

	RegionInfo* deriveRegion(const Region *region) const;
	RegionInfo* copyWithoutMesh() const;

	void prepare(const std::vector<Element*> &region, size_t begin, size_t end);

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

private:
	DistributedInfo(const Mesh *mesh, double domainShrinkRatio, double clusterShrinkRatio);
};

}
}



#endif /* SRC_OUTPUT_DISTRIBUTEDINFO_H_ */
