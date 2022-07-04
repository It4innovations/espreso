
#ifndef SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_

#include "input/input.h"

namespace espreso {

class InputConfiguration;

class InputOpenFoam: public Input {
public:
	InputOpenFoam();
	~InputOpenFoam();

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	void variables(Mesh &mesh);

protected:
	InputMesh<OrderedNodes, OrderedFaces, OrderedRegions> mesh;
	InputOpenFoam *loader;
};

class InputOpenFoamSequential: public InputOpenFoam {
public:
	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	void variables(Mesh &mesh);
};

class InputOpenFoamParallel: public InputOpenFoam {
public:
	InputOpenFoamParallel(esint domains): domains(domains) {}
	esint domains;

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	void variables(Mesh &mesh);
};

}

#endif /* SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_ */
