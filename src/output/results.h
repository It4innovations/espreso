
#ifndef SRC_OUTPUT_RESULTS_H_
#define SRC_OUTPUT_RESULTS_H_

#include "esmesh.h"

namespace espreso {
namespace output {

class Results {

public:
	virtual void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shringCluster) = 0;

	virtual ~Results() {};

protected:
	Results(const Mesh &mesh, const std::string &path): _mesh(mesh), _path(path) {};

	const Mesh &_mesh;
	std::string _path;
};

}
}



#endif /* SRC_OUTPUT_RESULTS_H_ */
