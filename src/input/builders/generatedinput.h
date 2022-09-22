
#ifndef SRC_INPUT_GENERATEDINPUT_H_
#define SRC_INPUT_GENERATEDINPUT_H_

#include <input/builders/inputold.h>

namespace espreso {

class GeneratedInput: public InputOLD {
public:
	GeneratedInput(MeshBuilder &dMesh, bool needSynchronization);

protected:
	void removeDanglingNodes();
	void synchronizeGlobalIndices();
};

}


#endif /* SRC_INPUT_GENERATEDINPUT_H_ */
