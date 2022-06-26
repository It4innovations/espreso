
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
	InputMesh<OrderedUniqueNodes, OrderedUniqueFaces, OrderedRegions> mesh;
};

}

#endif /* SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_ */
