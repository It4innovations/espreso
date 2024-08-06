
#include "coupling.h"
#include "config/configuration.hpp"

using namespace espreso;

CouplingSettings::CouplingSettings()
{
    REGISTER(configuration, ECFMetaData()
            .setdescription({ "PreCICE configuration file." })
            .setdatatype({ ECFDataType::STRING }));
    REGISTER(solver, ECFMetaData()
            .setdescription({ "Solver name." })
            .setdatatype({ ECFDataType::STRING }));
    REGISTER(mesh, ECFMetaData()
            .setdescription({ "Mesh name." })
            .setdatatype({ ECFDataType::STRING }));
    REGISTER(data_in, ECFMetaData()
            .setdescription({ "Data to be read from other solver." })
            .setdatatype({ ECFDataType::STRING }));
    REGISTER(data_out, ECFMetaData()
            .setdescription({ "Data to be write from espreso." })
            .setdatatype({ ECFDataType::STRING }));
}

CouplingConfiguration::CouplingConfiguration()
{
    active = false;
    REGISTER(active, ECFMetaData()
            .setdescription({ "Activate coupling with PreCICE." })
            .setdatatype({ ECFDataType::BOOL }));

    REGISTER(dummy, ECFMetaData()
            .setdescription({ "Settings for dummy-coupler." }));
}
