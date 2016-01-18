
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"

namespace esinput {

class OpenFOAM: public ExternalLoader {

public:
	OpenFOAM(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

	void open() {};
	void close() {};

private:
	std::string _path;
};

}



#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
