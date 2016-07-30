

#ifndef SRC_OUTPUT_PARAVIEW_PARAVIEW_H_
#define SRC_OUTPUT_PARAVIEW_PARAVIEW_H_

#include "../results.h"

namespace espreso {
namespace output {

class Paraview: public Results {

public:
	Paraview(const Mesh &mesh, const std::string &path);

	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);
};

}
}



#endif /* SRC_OUTPUT_PARAVIEW_PARAVIEW_H_ */
