
#ifndef SRC_OUTPUT_GENERIC_VTK_H_
#define SRC_OUTPUT_GENERIC_VTK_H_

#include "../store.h"

namespace espreso {
namespace output {

class Generic: public ResultStore {

public:
	Generic(const Mesh &mesh, const std::string &path);

	void store(double shrinkSubdomain, double shrinkCluster);
	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);
};

}
}

#endif /* SRC_OUTPUT_GENERIC_VTK_H_ */
