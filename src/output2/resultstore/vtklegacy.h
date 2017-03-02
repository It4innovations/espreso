
#ifndef SRC_OUTPUT2_RESULTSTORE_VTKLEGACY_H_
#define SRC_OUTPUT2_RESULTSTORE_VTKLEGACY_H_

#include "../resultstore.h"

namespace espreso {

namespace output {

class VTKLegacy: public ResultStore {

public:
	VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual ~VTKLegacy() {};

protected:
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data);
};

}
}


#endif /* SRC_OUTPUT2_RESULTSTORE_VTKLEGACY_H_ */
