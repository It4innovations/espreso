
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"
#include "foam/foamfile.h"

namespace esinput {

class OpenFOAM: public ExternalLoader {

public:
	OpenFOAM(const Options &options, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries);

	void open() {};
	void close() {};

private:

	ParseError* computePolyMeshPath(int rank, int size);
	void solveParseError(ParseError *error);

	std::string _projectPath;
	std::string _polyMeshPath;
	int _rank;
	int _size;
};

}



#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
