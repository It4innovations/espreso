
#ifndef INPUT_ANSYS_ANSYS_H_
#define INPUT_ANSYS_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>

#include "../loader.h"


namespace esinput {

class Ansys: public ExternalLoader {

public:
	Ansys(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

private:
	mesh::Element* createElement(eslocal *indices, eslocal n);

	static size_t getLinesCount(const std::string &file);
	struct TestEOL {
		bool operator()(char c) {
			return c == '\n';
		}
	};

	std::string _path;
};

}





#endif /* INPUT_ANSYS_ANSYS_H_ */
