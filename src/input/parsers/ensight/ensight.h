
#ifndef SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_
#define SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_

#include "input/input.h"

namespace espreso {

class InputConfiguration;
struct EnsightData;

class InputEnsight: public Input {
public:
	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	void variables(Mesh &mesh);

	~InputEnsight();

protected:
	EnsightData *data;
	InputMesh<OrderedNodes, OrderedElements, OrderedRegions> mesh;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_ */
