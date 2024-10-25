
#ifndef SRC_CONFIG_ECF_COUPLING_H_
#define SRC_CONFIG_ECF_COUPLING_H_

#include "config/description.h"

#include <string>

namespace espreso {

struct CouplingDataInConfiguration: public ECFDescription {
    bool force, pressure, stress;

    CouplingDataInConfiguration();

    bool isactive()
    {
        return force || pressure || stress;
    }
};

struct CouplingDataOutConfiguration: public ECFDescription {
    bool displacement, velocity;

    CouplingDataOutConfiguration();

    bool isactive()
    {
        return displacement || velocity;
    }
};

struct CouplingSettings: public ECFDescription {
    std::string configuration, solver, mesh;

    CouplingDataInConfiguration data_in;
    CouplingDataOutConfiguration data_out;

    CouplingSettings();
};

struct CouplingConfiguration: public CouplingSettings {
    CouplingSettings dummy;

    CouplingConfiguration();

    bool isactive()
    {
        return data_in.isactive() || data_out.isactive();
    }
};

}

#endif /* SRC_CONFIG_ECF_COUPLING_H_ */
