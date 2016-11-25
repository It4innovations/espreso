
#ifndef INPUT_ANSYS_ANSYS_H_
#define INPUT_ANSYS_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>


#include "../loader.h"
#include "utils.h"
#include "parser.h"

namespace espreso {
namespace input {

class AnsysWorkbench: public Loader {

public:
	static void load(Mesh &mesh, const ArgsConfiguration &configuration, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from Ansys/Workbench format from file " << configuration.path;
		AnsysWorkbench workbench(mesh, configuration, rank, size);
		workbench.fill();
	}

protected:
	AnsysWorkbench(Mesh &mesh, const ArgsConfiguration &configuration, int rank, int size)
	: Loader(mesh), _path(configuration.path), _parser(mesh) { };

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void materials(std::vector<Material> &materials);
	void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours);

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





#endif /* INPUT_ANSYS_ANSYS_H_ */
