
#ifndef SRC_INPUT_GENERATEDINPUT_H_
#define SRC_INPUT_GENERATEDINPUT_H_

#include "input.h"

namespace espreso {

class GeneratedInput: public Input {
public:
	GeneratedInput(MeshBuilder &dMesh, bool needSynchronization);

protected:
	void removeDanglingNodes();
	void synchronizeGlobalIndices();
};

}


#endif /* SRC_INPUT_GENERATEDINPUT_H_ */
