
#include "transientfirstorderimplicit.h"
#include "config/configuration.hpp"

espreso::AutoTimeSteppingConfiguration::AutoTimeSteppingConfiguration()
{
    allowed = false;
    REGISTER(allowed, ECFMetaData()
            .setdescription({ "Allow" })
            .setdatatype({ ECFDataType::BOOL }));

    ecfdescription->addSeparator();

    min_time_step = 1e-3;
    max_time_step = 1;
    REGISTER(min_time_step, ECFMetaData()
            .setdescription({ "Minimum time step" })
            .setdatatype({ ECFDataType::FLOAT }));
    REGISTER(max_time_step, ECFMetaData()
            .setdescription({ "Maximum time step" })
            .setdatatype({ ECFDataType::FLOAT }));

    ecfdescription->addSeparator();

    oscilation_limit = 0.5;
    IDFactor = 3;
    REGISTER(oscilation_limit, ECFMetaData()
            .setdescription({ "Oscilation limit" })
            .setdatatype({ ECFDataType::FLOAT }));
    REGISTER(IDFactor, ECFMetaData()
            .setdescription({ "I/D factor" })
            .setdatatype({ ECFDataType::FLOAT }));

    points_per_period = 0;
    REGISTER(points_per_period, ECFMetaData()
            .setdescription({ "Minimal points per period" })
            .setdatatype({ ECFDataType::INTEGER }));
}

espreso::TransientFirstOrderImplicitSolverConfiguration::TransientFirstOrderImplicitSolverConfiguration()
{
    method = METHOD::CRANK_NICOLSON;
    REGISTER(method, ECFMetaData()
            .setdescription({ "Method" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("CRANK_NICOLSON").setdescription("Alpha = 0.5."))
            .addoption(ECFOption().setname("FORWARD_DIFF").setdescription("Alpha = ??."))
            .addoption(ECFOption().setname("GALERKIN").setdescription("Alpha = 2 / 3."))
            .addoption(ECFOption().setname("BACKWARD_DIFF").setdescription("Alpha = 1."))
            .addoption(ECFOption().setname("USER").setdescription("User defined Alpha from interval <0, 1).")));

    alpha = 0.5;
    REGISTER(alpha, ECFMetaData()
            .setdescription({ "Alpha" })
            .setdatatype({ ECFDataType::FLOAT })
            .allowonly([&] () { return method == METHOD::USER; }));

    ecfdescription->addSpace()->metadata.noexport();

    time_step = 0.1;
    REGISTER(time_step, ECFMetaData()
            .setdescription({ "Time step" })
            .setdatatype({ ECFDataType::FLOAT }));

    REGISTER(auto_time_stepping, ECFMetaData()
            .setdescription({ "Auto time stepping" })
            .setcollapsed());


}



