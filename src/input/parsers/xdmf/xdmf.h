
#ifndef SRC_INPUT_FORMATS_XDMF_XDMF_H_
#define SRC_INPUT_FORMATS_XDMF_XDMF_H_

#include "input/meshbuilder.h"

namespace espreso {

struct InputConfiguration;

class XDMFLoader: public MeshBuilder {
public:
    XDMFLoader(const InputConfiguration &configuration);
    void load();

protected:
    const InputConfiguration &_configuration;
};

}



#endif /* SRC_INPUT_FORMATS_XDMF_XDMF_H_ */
