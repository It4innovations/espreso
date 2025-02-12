
#include "nonlinear.h"
#include "config/configuration.hpp"

espreso::NonLinearSolverConfiguration::NonLinearSolverConfiguration(const std::string &firstResidualName, const std::string &secondResidualName)
{
    method = METHOD::NEWTON_RAPHSON;
    REGISTER(method, ECFMetaData()
            .setdescription({ "Method" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("NEWTON_RAPHSON").setdescription("Newton-Raphson."))
            .addoption(ECFOption().setname("MODIFIED_NEWTON_RAPHSON").setdescription("Newton-Raphson without re-assembling of stiffness matrices.")));

    ecfdescription->addSpace()->metadata.noexport();

    check_first_residual = true;
    check_second_residual = false;
    ecfdescription->registerParameter("check_" + firstResidualName, check_first_residual, ECFMetaData()
            .setdescription({ "Temperature convergance" })
            .setdatatype({ ECFDataType::FLOAT }));
    ecfdescription->registerParameter("check_" + secondResidualName, check_second_residual, ECFMetaData()
            .setdescription({ "Heat convergance" })
            .setdatatype({ ECFDataType::FLOAT }));

    requested_first_residual = requested_second_residual = 1e-3;
    ecfdescription->registerParameter("requested_" + firstResidualName + "_residual", requested_first_residual, ECFMetaData()
            .setdescription({ "Temperature precision" })
            .setdatatype({ ECFDataType::FLOAT }));
    ecfdescription->registerParameter("requested_" + secondResidualName + "_residual", requested_second_residual, ECFMetaData()
            .setdescription({ "Heat precision" })
            .setdatatype({ ECFDataType::FLOAT }));

    ecfdescription->addSeparator();

    stepping = STEPPINGG::FALSE;
    REGISTER(stepping, ECFMetaData()
            .setdescription({ "Substepping" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("TRUE").setdescription("Turn on."))
            .addoption(ECFOption().setname("FALSE").setdescription("Turn off."))
            .addoption(ECFOption().setname("AUTO").setdescription("Automatic.")));

    substeps = 1;
    REGISTER(substeps, ECFMetaData()
            .setdescription({ "Number of substeps" })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    max_iterations = 15;
    REGISTER(max_iterations, ECFMetaData()
            .setdescription({ "Maximum iterations" })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    ecfdescription->addSpace()->metadata.noexport();

    line_search = tangent_matrix_correction = adaptive_precision = false;
    REGISTER(line_search, ECFMetaData()
            .setdescription({ "Line search" })
            .setdatatype({ ECFDataType::BOOL }));
    REGISTER(tangent_matrix_correction, ECFMetaData()
            .setdescription({ "Tangens stiffness matrix" })
            .setdatatype({ ECFDataType::BOOL }));
    REGISTER(adaptive_precision, ECFMetaData()
            .setdescription({ "Adaptive precision" })
            .setdatatype({ ECFDataType::BOOL }));

    ecfdescription->addSpace()->metadata.noexport();

    r_tol = 0.1;
    c_fact = 0.8;
    REGISTER(r_tol, ECFMetaData()
            .setdescription({ "R-tolerance" })
            .setdatatype({ ECFDataType::FLOAT }));
    REGISTER(c_fact, ECFMetaData()
            .setdescription({ "C-factor" })
            .setdatatype({ ECFDataType::FLOAT }));
}
