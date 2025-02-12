
#ifndef SRC_INPUT_FORMATS_VTKLEGACY_VTKLEGACY_H_
#define SRC_INPUT_FORMATS_VTKLEGACY_VTKLEGACY_H_

#include "input/meshbuilder.h"

namespace espreso {

struct InputConfiguration;

class VTKLegacyLoader: public MeshBuilder {
public:
    VTKLegacyLoader(const InputConfiguration &configuration);
    void load();

protected:
    const InputConfiguration &_configuration;
};

}

#endif /* SRC_INPUT_FORMATS_VTKLEGACY_VTKLEGACY_H_ */
