
#ifndef SRC_WRAPPERS_CATALYST_W_CATALYST_H_
#define SRC_WRAPPERS_CATALYST_W_CATALYST_H_

class vtkUnstructuredGrid;
class vtkCPProcessor;
class vtkCPDataDescription;

namespace espreso {

class Mesh;

class Catalyst {

public:
	static constexpr bool islinked();

	Catalyst();
	~Catalyst();

	void update();

protected:
	vtkCPProcessor *_processor;
	vtkUnstructuredGrid *_VTKGrid;
	vtkCPDataDescription *_dataDescription;
};

}

#endif /* SRC_WRAPPERS_CATALYST_W_CATALYST_H_ */
