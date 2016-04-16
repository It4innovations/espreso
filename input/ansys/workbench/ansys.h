
#ifndef INPUT_ANSYS_WORKBENCH_ANSYS_H_
#define INPUT_ANSYS_WORKBENCH_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "../../loader.h"
#include "../utils.h"

namespace espreso {
namespace input {

class AnsysWorkbench: public Loader {

public:
	static void load(Mesh &mesh, const Options &options, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from Ansys/Workbench format from file " << options.path;
		AnsysWorkbench workbench(mesh, options, rank, size);
		workbench.fill();
	}

protected:
	AnsysWorkbench(Mesh &mesh, const Options &options, int rank, int size)
	: Loader(mesh), _path(options.path) { };

	void points(Coordinates &coordinates, size_t &DOFs);
	void elements(std::vector<Element*> &elements);
	void boundaryConditions(Coordinates &coordinates);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

	void open()
	{
		_file.open(_path.c_str());
		if (!_file.is_open()) {
			ESINFO(ERROR) << "Cannot load mesh from file: " << _path;
		}
	}

	void close()
	{
		_file.close();
	}

private:
	std::ifstream _file;
	std::string _path;
};

}
}





#endif /* INPUT_ANSYS_WORKBENCH_ANSYS_H_ */
