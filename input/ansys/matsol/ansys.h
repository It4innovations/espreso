
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
	AnsysMatsol(const Options &options, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void faces(mesh::Faces &faces) {};
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries, std::vector<int> &neighbours);

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
