

#ifndef SRC_OUTPUT_PARAVIEW_PARAVIEW_H_
#define SRC_OUTPUT_PARAVIEW_PARAVIEW_H_

#include "../store.h"

namespace espreso {
namespace output {

class Paraview: public ResultStore {

public:
	Paraview(const Mesh &mesh, const std::string &path);

	void store(double shrinkSubdomain, double shrinkCluster);
	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);
};

}
}



#endif /* SRC_OUTPUT_PARAVIEW_PARAVIEW_H_ */
