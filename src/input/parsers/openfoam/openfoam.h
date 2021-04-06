
#ifndef SRC_INPUT_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_OPENFOAM_OPENFOAM_H_

#include "input/meshbuilder.h"
#include "basis/io/inputfile.h"

#include <cstddef>
#include <string>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;
struct OpenFOAMSet;

struct OpenFOAMData: public MeshBuilder {
	esint nelements;
	std::vector<esint> fIDs, fsize, fnodes, owner, neighbor;
};

class OpenFOAMLoader: public OpenFOAMData {

public:
	static std::string cellprefix;
	OpenFOAMLoader(const InputConfiguration &configuration);
	void load();

protected:
	void readData();
	void parseData();

	void buildFaces();

	void collectFaces();
	void buildElements();

	const InputConfiguration &_configuration;

	InputFile _points, _faces, _owner, _neighbor, _boundary;
	InputFile _pointZones, _faceZones, _cellZones;

	std::vector<OpenFOAMSet> _sets;
	std::vector<esint> _fdist, _edist;
};

}



#endif /* SRC_INPUT_OPENFOAM_OPENFOAM_H_ */
