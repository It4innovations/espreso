
#ifndef SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_
#define SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_

#include "../resultstore.h"

namespace espreso {

namespace output {

class VTKLegacy: public ResultStore {

public:
	VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual ~VTKLegacy() {};

protected:
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const DataArrays &data);
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution);

	// Legacy format cannot be linked
	virtual void linkClusters(const std::string &root, const std::string &name, const DataArrays &data) {};
	virtual void linkClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution) {};

	virtual void linkSteps(const std::string &name, const std::vector<std::pair<std::string, Step> > &steps) {};
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_ */
