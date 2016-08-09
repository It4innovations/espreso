
#ifndef SRC_OUTPUT_VTK_VTK_H_
#define SRC_OUTPUT_VTK_VTK_H_

#include "../results.h"
#include "../../assembler/constraints/equalityconstraints.h"

namespace espreso {
namespace output {

class VTK: public Results {

public:
	VTK(const Mesh &mesh, const std::string &path);

	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);

	static void gluing(const Mesh &mesh, const EqualityConstraints &constraints, const std::string &path, double shrinkSubdomain, double shrinkCluster);

	static void mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster);
	static void properties(const Mesh &mesh, const std::string &path, std::vector<Property> properties, double shrinkSubdomain, double shrinkCluster);

	static void fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
};

}
}

#endif /* SRC_OUTPUT_VTK_VTK_H_ */
