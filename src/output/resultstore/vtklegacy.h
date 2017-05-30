
#ifndef SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_
#define SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_

#include "../resultstore.h"

namespace espreso {
namespace output {

class VTKLegacy: public ResultStore {

public:
	VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path, MeshInfo::InfoMode mode = MeshInfo::EMPTY);
	virtual ~VTKLegacy() {};

protected:
	virtual std::string store(const std::string &name, const RegionData &regionData);

	// Legacy format cannot be linked
	virtual std::string linkClusters(const std::string &root, const std::string &name, const RegionData &regionData) { return ""; };
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps) {};
};

}
}


#endif /* SRC_OUTPUT_RESULTSTORE_VTKLEGACY_H_ */
