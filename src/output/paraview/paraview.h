

#ifndef SRC_OUTPUT_PARAVIEW_PARAVIEW_H_
#define SRC_OUTPUT_PARAVIEW_PARAVIEW_H_

#include "../store.h"

namespace espreso {
namespace output {

class Paraview: public Store {

public:
	Paraview(const Mesh &mesh, const std::string &path);

	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType) { ESINFO(GLOBAL_ERROR) << "Implement store property"; }
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType) { ESINFO(GLOBAL_ERROR) << "Implement store property"; }

	void store(std::vector<std::vector<double> > &displacement, double shrinkSubdomain, double shrinkCluster);
};

}
}



#endif /* SRC_OUTPUT_PARAVIEW_PARAVIEW_H_ */
