
#ifndef SRC_CONFIG_ECF_COUPLING_H_
#define SRC_CONFIG_ECF_COUPLING_H_

#include "config/description.h"

#include <string>
#include <map>

namespace espreso {

struct DataExchangeConfiguration: public ECFDescription {
    struct DataRead: public ECFDescription {
        bool force, pressure, stress;

        DataRead(): force(false), pressure(false), stress(false) {}
    };

    struct DataWrite: public ECFDescription {
        bool displacement, velocity;

        DataWrite(): displacement(false), velocity(false) {}
    };

    bool direct;
    bool centers;
    DataRead read;
    DataWrite write;

    DataExchangeConfiguration();

    bool isactive() const
    {
        return read.force || read.pressure || read.stress || write.displacement || write.velocity;
    }
};

struct CouplingConfiguration: public ECFDescription {

    std::string configuration, solver, mesh, centers;

    std::map<std::string, DataExchangeConfiguration> exchange;

    CouplingConfiguration();

    bool isactive() const
    {
        for (auto data = exchange.cbegin(); data != exchange.cend(); ++data) {
            if (data->second.isactive()) {
                return true;
            }
        }
        return false;
    }
};

}

#endif /* SRC_CONFIG_ECF_COUPLING_H_ */
