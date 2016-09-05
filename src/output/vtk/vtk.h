
#ifndef SRC_OUTPUT_VTK_VTK_H_
#define SRC_OUTPUT_VTK_VTK_H_

#include "../../assembler/constraints/equalityconstraints.h"
#include "../store.h"

namespace espreso {
namespace output {

class VTK: public Store {

public:
	VTK(const Mesh &mesh, const std::string &path, double shrinkSubdomain = config::output::SUBDOMAINS_SHRINK_RATIO, double shringCluster = config::output::CLUSTERS_SHRINK_RATIO);

	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType);
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);

	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);

	static void gluing(const Mesh &mesh, const EqualityConstraints &constraints, const std::string &path, size_t dofs, double shrinkSubdomain, double shrinkCluster);

	static void mesh(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shrinkCluster);
	static void properties(const Mesh &mesh, const std::string &path, std::vector<Property> properties, double shrinkSubdomain, double shrinkCluster);

	static void fixPoints(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
	static void corners(const Mesh &mesh, const std::string &path, double shrinkSubdomain, double shringCluster);
};

}
}

#endif /* SRC_OUTPUT_VTK_VTK_H_ */
