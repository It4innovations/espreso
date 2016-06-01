
#ifndef INPUT_ANSYS_MATSOL_ANSYS_H_
#define INPUT_ANSYS_MATSOL_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>

#include "../../loader.h"
#include "../utils.h"

namespace espreso {
namespace input {

class AnsysMatsol: public Loader {

public:
	static void load(Mesh &mesh, const Configuration &configuration, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from Ansys/Matsol format from directory " << configuration.path;

		AnsysMatsol matsol(mesh, configuration, rank, size);
		matsol.fill();
	}

protected:
	AnsysMatsol(Mesh &mesh, const Configuration &configuration, int rank, int size)
	: Loader(mesh), _path(configuration.path) { };

	void points(Coordinates &coordinates, size_t &DOFs);
	void elements(std::vector<Element*> &elements);
	void materials(std::vector<Material> &materials);
	void boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions);
	void clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours);

private:
	static size_t getLinesCount(const std::string &file);
	struct TestEOL {
		bool operator()(char c) {
			return c == '\n';
		}
	};

	std::string _path;
};

}
}




#endif /* INPUT_ANSYS_MATSOL_ANSYS_H_ */
