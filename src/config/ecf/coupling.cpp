
#include "coupling.h"
#include "config/configuration.hpp"

using namespace espreso;

DataExchangeConfiguration::DataExchangeConfiguration()
{
    direct = false;
    REGISTER(direct, ECFMetaData()
            .setdescription({ "Direct read of data." })
            .setdatatype({ ECFDataType::BOOL }));

    centers = false;
    REGISTER(centers, ECFMetaData()
            .setdescription({ "Exchange data on edge/face centers." })
            .setdatatype({ ECFDataType::BOOL }));

    read.force = false;
    ecfdescription->registerParameter("force", read.force, ECFMetaData()
            .setdescription({ "Read FORCE." })
            .setdatatype({ ECFDataType::BOOL }));
    read.pressure = false;
    ecfdescription->registerParameter("pressure", read.pressure, ECFMetaData()
            .setdescription({ "Read PRESSURE." })
            .setdatatype({ ECFDataType::BOOL }));
    read.stress = false;
    ecfdescription->registerParameter("stress", read.stress, ECFMetaData()
            .setdescription({ "Read STRESS." })
            .setdatatype({ ECFDataType::BOOL }));

    write.displacement = false;
    ecfdescription->registerParameter("displacement", write.displacement, ECFMetaData()
            .setdescription({ "Write DISPLACEMENT." })
            .setdatatype({ ECFDataType::BOOL }));
    write.velocity = false;
    ecfdescription->registerParameter("velocity", write.velocity, ECFMetaData()
            .setdescription({ "Write VELOCITY." })
            .setdatatype({ ECFDataType::BOOL }));
}

CouplingConfiguration::CouplingConfiguration()
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
    REGISTER(centers, ECFMetaData()
            .setdescription({ "Mesh centers name." })
            .setdatatype({ ECFDataType::STRING }));

    REGISTER(exchange, ECFMetaData()
            .setdescription({ "Mesh name.", "Data settings" })
            .setdatatype({ ECFDataType::STRING })
            .setpattern({ "Mesh-Solid" }));
}
