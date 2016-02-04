
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

	/** @brief Project path. */
	std::string _projectPath;

	/** @brief Path to PolyMesh, it contains also rank number for divided cases. */
	std::string _polyMeshPath;

	/** @brief Assigned rank, 0 for non MPI runs.*/
	int _rank;

	/** @brief Number of processes, 1 for non MPI runs*/
	int _size;
};

}



#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
