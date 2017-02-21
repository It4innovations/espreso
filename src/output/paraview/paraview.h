

#ifndef SRC_OUTPUT_PARAVIEW_PARAVIEW_H_
#define SRC_OUTPUT_PARAVIEW_PARAVIEW_H_

#include "../resultstore.h"

namespace espreso {
namespace store {

class Paraview: public ResultStore {

public:
	Paraview(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);

	virtual void storeGeometry(size_t timeStep = -1);
	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType);
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);
	virtual void finalize();
};

}
}



#endif /* SRC_OUTPUT_PARAVIEW_PARAVIEW_H_ */
