
#ifndef INPUT_ANSYS_MATSOL_ANSYS_H_
#define INPUT_ANSYS_MATSOL_ANSYS_H_

#include <string>
#include <sstream>
#include <fstream>

#include "../../loader.h"
#include "../utils.h"


namespace esinput {

class AnsysMatsol: public ExternalLoader {

public:
	AnsysMatsol(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

	void open() {};
	void close() {};

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





#endif /* INPUT_ANSYS_MATSOL_ANSYS_H_ */
