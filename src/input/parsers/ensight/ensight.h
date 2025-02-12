
#ifndef SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_
#define SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_

#include "input/meshbuilder.h"

namespace espreso {

struct InputConfiguration;

class EnsightLoader: public MeshBuilder {
public:
    EnsightLoader(const InputConfiguration &configuration);
    void load();

protected:
    const InputConfiguration &_configuration;
};

}

#endif /* SRC_INPUT_FORMATS_ENSIGHT_ENSIGHT_H_ */
