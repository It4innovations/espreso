
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"
#include "foam/foamfile.h"
#include "foam/face.h"
#include "foam/dictionary.h"
#include "foam/elementbuilder.h"
#include "foam/cellzone.h"

namespace esinput {

class OpenFOAM: public ExternalLoader {

public:
	OpenFOAM(const Options &options, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);
	void faces(mesh::Faces &faces);
	void boundaryConditions(mesh::Coordinates &coordinates);
	void clusterBoundaries(mesh::Mesh &mesh, mesh::Boundaries &boundaries, std::vector<int> &neighbours);

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

	/** @brief Temporary storage for faces*/
	std::vector<Face> _faces;

	/** @brief Temporary storage for cell zones*/
	std::vector<CellZone> _cellZones;
};

}

#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
