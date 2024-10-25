
#include "coupling.h"
#include "config/configuration.hpp"

using namespace espreso;

CouplingDataInConfiguration::CouplingDataInConfiguration()
{
    force = false;
    REGISTER(force, ECFMetaData()
            .setdescription({ "Read FORCE." })
            .setdatatype({ ECFDataType::BOOL }));
    pressure = false;
    REGISTER(pressure, ECFMetaData()
            .setdescription({ "Read PRESSURE." })
            .setdatatype({ ECFDataType::BOOL }));
    stress = false;
    REGISTER(stress, ECFMetaData()
            .setdescription({ "Read STRESS." })
            .setdatatype({ ECFDataType::BOOL }));
}

CouplingDataOutConfiguration::CouplingDataOutConfiguration()
{
    displacement = false;
    REGISTER(displacement, ECFMetaData()
            .setdescription({ "Write DISPLACEMENT." })
            .setdatatype({ ECFDataType::BOOL }));
    velocity = false;
    REGISTER(velocity, ECFMetaData()
            .setdescription({ "Write VELOCITY." })
            .setdatatype({ ECFDataType::BOOL }));
}

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
            .setdescription({ "Read data." }));
    REGISTER(data_out, ECFMetaData()
            .setdescription({ "Writted data." }));
}

CouplingConfiguration::CouplingConfiguration()
{
    REGISTER(dummy, ECFMetaData()
            .setdescription({ "Settings for dummy-coupler." }));
}
