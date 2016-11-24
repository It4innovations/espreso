
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

	virtual ~ClusterGenerator() {};

protected:
	ClusterGenerator(Mesh &mesh): _mesh(mesh) {};

	Mesh &_mesh;
};

}
}



#endif /* SRC_INPUT_GENERATOR_GENERATOR_H_ */
