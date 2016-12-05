
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"
#include "foam/foamfile.h"
#include "foam/face.h"
#include "foam/dictionary.h"
#include "foam/elementbuilder.h"
#include "foam/cellzone.h"
#include "foam/Boundary.h"

#include "../../config/description.h"

namespace espreso {
namespace input {

class OpenFOAM: public Loader {

public:
	static void load(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size)
	{
		ESINFO(OVERVIEW) << "Load mesh from OpenFOAM format from directory " << configuration.openfoam.path;
		OpenFOAM openfoam(configuration, mesh, rank, size);
		openfoam.fill();
	}

	bool faceBased() const { return true; }

protected:
	OpenFOAM(const GlobalConfiguration &configuration, Mesh &mesh, int rank, int size);

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
		mesh.partitiate(_openfoam.domains);
		return true;
	}

private:

	ParseError* computePolyMeshPath(int rank, int size);
	void solveParseError(ParseError *error);

	ESPRESOInput _openfoam;

	/** @brief Project path. */
	std::string _projectPath;

	/** @brief Path to PolyMesh, it contains also rank number for divided cases. */
	std::string _polyMeshPath;

	/** @brief Assigned rank, 0 for non MPI runs.*/
	int _rank;

	/** @brief Number of processes, 1 for non MPI runs*/
	int _size;

	/** @brief Temporary storage for loaded boundaries*/
	std::vector<Boundary > _boundaries;
};

}
}

#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
