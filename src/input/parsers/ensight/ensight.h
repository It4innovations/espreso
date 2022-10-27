
#ifndef SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_
#define SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_

#include "input/input.h"

namespace espreso {

class InputConfiguration;
struct EnsightData;

class InputEnsight: public Input {
public:
	InputEnsight();

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh);
	void finishVariables();
	int timeSteps();
	void nextTimeStep(Mesh &mesh);

	~InputEnsight();

protected:
	EnsightData *data;

	NodesBlocks nodes;
	ElementsBlocks elements;
	VariablesBlocks variables;
	OrderedRegions regions;

	int timeset;
	size_t timestep;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_ */
