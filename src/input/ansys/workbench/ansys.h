
#ifndef INPUT_ANSYS_WORKBENCH_ANSYS_H_
#define INPUT_ANSYS_WORKBENCH_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>


#include "../../loader.h"
#include "../utils.h"
#include "parser.h"

namespace espreso {
namespace input {

class AnsysWorkbench: public Loader {

public:
	static void load(Mesh &mesh, const Configuration &configuration, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from Ansys/Workbench format from file " << configuration.path;
		AnsysWorkbench workbench(mesh, configuration, rank, size);
		workbench.fill();
	}

protected:
	AnsysWorkbench(Mesh &mesh, const Configuration &configuration, int rank, int size)
	: Loader(mesh), _path(configuration.path) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials);
	void boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

	void open()
	{
		_parser.open(_path);
	}

	void close()
	{
		_parser.close();
	}

private:
	std::string _path;
	WorkbenchParser _parser;
};

}
}





#endif /* INPUT_ANSYS_WORKBENCH_ANSYS_H_ */
