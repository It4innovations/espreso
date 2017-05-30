
#ifndef SRC_OUTPUT_RESULTSTORE_CATALYST_H_
#define SRC_OUTPUT_RESULTSTORE_CATALYST_H_

#include "vtkxmlascii.h"

class vtkCPProcessor;
class vtkCPDataDescription;
class vtkFieldData;

namespace espreso {

namespace output {

class Catalyst: public VTKXMLASCII {

public:
	Catalyst(const OutputConfiguration &output, const Mesh *mesh, const std::string &path, MeshInfo::InfoMode mode = MeshInfo::PREPARE);
	virtual ~Catalyst();

	virtual void storeSettings(const Step &step);
	virtual void storeSettings(size_t steps);
	virtual void storeSettings(const std::vector<size_t> &steps);
	virtual void storeFETIData(const Step &step, const Instance &instance);

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution);

protected:
	vtkCPProcessor *_processor;
	vtkCPDataDescription *_dataDescription;
	vtkFieldData *_fieldData;
};

}
}



#endif /* SRC_OUTPUT_RESULTSTORE_CATALYST_H_ */
