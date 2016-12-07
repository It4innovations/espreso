
#ifndef INPUT_ANSYS_ANSYS_H_
#define INPUT_ANSYS_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "../loader.h"
#include "utils.h"
#include "parser.h"

#include "../../config/description.h"

namespace espreso {
namespace input {

class AnsysWorkbench: public Loader {

public:
	static void load(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from Ansys/Workbench format from file " << configuration.workbench.path;
		AnsysWorkbench workbench(configuration, mesh, rank, size);
		workbench.fill();
	}

protected:
	AnsysWorkbench(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size)
	: Loader(configuration, mesh), _workbench(configuration.workbench), _parser(mesh) { };

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
	bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners)
	{
		mesh.partitiate(_workbench.domains);
		return true;
	}

	void open()
	{
		_parser.open(_workbench.path);
	}

	void close()
	{
		_parser.close();
	}

private:
	ESPRESOInput _workbench;
	WorkbenchParser _parser;
};

}
}





#endif /* INPUT_ANSYS_ANSYS_H_ */
