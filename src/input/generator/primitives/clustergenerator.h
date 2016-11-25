
#ifndef SRC_INPUT_GENERATOR_GENERATOR_H_
#define SRC_INPUT_GENERATOR_GENERATOR_H_

#include "../../meshgenerator/elements/elements.h"

namespace espreso {
namespace input {

class ClusterGenerator {

public:
	virtual void points(std::vector<Point> &points) =0;
	virtual void elements(std::vector<Element*> &elements) =0;
	virtual void boundaries(std::vector<Element*> &nodes, const std::vector<int> &neighbours) =0;

	virtual void uniformPartition(std::vector<eslocal> &partsPtrs, size_t subdomains) =0;
	virtual void uniformFixPoints(const std::vector<Element*> &nodes, std::vector<std::vector<Element*> > &fixPoints) =0;
	virtual void uniformCorners(const std::vector<Element*> &nodes, std::vector<Element*> &corners, size_t number, bool point, bool edge, bool face) =0;

	virtual ~ClusterGenerator() {};
};

}
}



#endif /* SRC_INPUT_GENERATOR_GENERATOR_H_ */
