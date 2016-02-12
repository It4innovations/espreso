
#ifndef INPUT_ANSYS_WORKBENCH_ANSYS_H_
#define INPUT_ANSYS_WORKBENCH_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "../../loader.h"
#include "../utils.h"


namespace esinput {

class AnsysWorkbench: public ExternalLoader {

public:
	AnsysWorkbench(const Options &options, size_t index, size_t size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

	void open()
	{
		_file.open(_path.c_str());
		if (!_file.is_open()) {
			ESLOG(eslog::ERROR) << "Cannot load mesh from file: " << _path;
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





#endif /* INPUT_ANSYS_WORKBENCH_ANSYS_H_ */
