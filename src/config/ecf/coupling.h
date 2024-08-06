
#ifndef SRC_CONFIG_ECF_COUPLING_H_
#define SRC_CONFIG_ECF_COUPLING_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct CouplingSettings: public ECFDescription {
    std::string configuration, solver, mesh, data_in, data_out;

    CouplingSettings();
};

struct CouplingConfiguration: public CouplingSettings {
    bool active;

    CouplingSettings dummy;

    CouplingConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_COUPLING_H_ */
