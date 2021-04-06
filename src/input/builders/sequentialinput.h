
#ifndef SRC_INPUT_SEQUENTIALINPUT_H_
#define SRC_INPUT_SEQUENTIALINPUT_H_

#include "input.h"

namespace espreso {

class SequentialInput: public Input {
public:
	SequentialInput(MeshBuilder &mesh);

protected:
	void removeDanglingNodes();
	void reindexNRegions();
	void reindexERegions();

	std::vector<esint> _dangling;
};

}


#endif /* SRC_INPUT_SEQUENTIALINPUT_H_ */
