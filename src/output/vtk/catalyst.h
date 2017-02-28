

#ifndef SRC_OUTPUT_VTK_CATALYST_H_
#define SRC_OUTPUT_VTK_CATALYST_H_

#include "vtk.h"

class vtkCPProcessor;
class vtkCPDataDescription;
class vtkFieldData;

namespace espreso {
namespace store {

class Catalyst: public VTK {

public:
	Catalyst(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);
	~Catalyst();

	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);
	virtual void finalize();

protected:
	vtkCPProcessor *processor;
	vtkCPDataDescription *dataDescription;
	vtkFieldData *fieldData;
	size_t timeStep;
};

}
}



#endif /* SRC_OUTPUT_VTK_CATALYST_H_ */
