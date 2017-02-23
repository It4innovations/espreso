

#ifndef SRC_OUTPUT_VTK_CATALYST_H_
#define SRC_OUTPUT_VTK_CATALYST_H_

#include "../resultstore.h"

class vtkUnstructuredGrid;
class vtkCPProcessor;
class vtkCPDataDescription;

namespace espreso {
namespace store {

class Catalyst: public ResultStore {

public:
	Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);
	~Catalyst();

	virtual void storeGeometry(size_t timeStep = -1);
	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType);
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);
	virtual void finalize();

protected:
	vtkUnstructuredGrid *VTKGrid;
	vtkCPProcessor *processor;
	vtkCPDataDescription *dataDescription;
	std::vector<void*> VTKDataArrays;
};

}
}



#endif /* SRC_OUTPUT_VTK_CATALYST_H_ */
