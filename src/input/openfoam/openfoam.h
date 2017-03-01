
#ifndef INPUT_OPENFOAM_OPENFOAM_H_
#define INPUT_OPENFOAM_OPENFOAM_H_

#include "../loader.h"
#include "foam/boundary.h"
#include "foam/foamfile.h"
#include "foam/face.h"
#include "foam/dictionary.h"
#include "foam/elementbuilder.h"
#include "foam/zones.h"


namespace espreso {

struct ESPRESOInput;

namespace input {

class OpenFOAM: public Loader {

public:
	static void load(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size);
	bool faceBased() const { return true; }

protected:
	OpenFOAM(const ESPRESOInput &configuration, Mesh &mesh, int rank, int size);

	void points(Coordinates &coordinates);
	void elements(std::vector<size_t> &bodies, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void materials(std::vector<Material> &materials) {};
	void regions(
				std::vector<Evaluator*> &evaluators,
				std::vector<Region*> &regions,
				std::vector<Element*> &elements,
				std::vector<Element*> &faces,
				std::vector<Element*> &edges,
				std::vector<Element*> &nodes);
	void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges);
	bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

private:

	ParseError* computePolyMeshPath(int rank, int size);
	void solveParseError(ParseError *error);

	const ESPRESOInput &_openfoam;

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
}

#endif /* INPUT_OPENFOAM_OPENFOAM_H_ */
