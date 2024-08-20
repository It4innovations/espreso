
#include "transientsecondorderimplicit.h"
#include "config/configuration.hpp"

using namespace espreso;

TransientSecondOrderImplicitSolverConfiguration::TransientSecondOrderImplicitSolverConfiguration()
{
	numerical_damping = 0;
	REGISTER(numerical_damping, ECFMetaData()
			.setdescription({ "Numerical Damping" })
			.setdatatype({ ECFDataType::FLOAT }));

	method = METHOD::NEWMARK;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Method" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NEWMARK").setdescription("NEWMARK method.")));

	alpha = 0.25;
	REGISTER(alpha, ECFMetaData()
			.setdescription({ "Alpha" })
			.setdatatype({ ECFDataType::FLOAT }));

	delta = 0.5;
	REGISTER(delta, ECFMetaData()
			.setdescription({ "Delta" })
			.setdatatype({ ECFDataType::FLOAT }));
	ecfdescription->addSpace();

	alphaM = 0;
    REGISTER(alphaM, ECFMetaData()
            .setdescription({ "alphaM" })
            .setdatatype({ ECFDataType::FLOAT }));
    ecfdescription->addSpace();

    alphaF = 0;
    REGISTER(alphaF, ECFMetaData()
            .setdescription({ "alphaF" })
            .setdatatype({ ECFDataType::FLOAT }));
    ecfdescription->addSpace();

	mass_matrix_type = MASS_MATRIX_TYPE::CONSISTENT;
	REGISTER(method, ECFMetaData()
			.setdescription({ "Method" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CONSISTENT").setdescription("Consistent method."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Diagonal method."))
			.addoption(ECFOption().setname("HRZDIAGONAL").setdescription("HRZ Diagonal method.")));

	ecfdescription->addSpace();

	time_step = 0.1;
	REGISTER(time_step, ECFMetaData()
			.setdescription({ "Time step" })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(auto_time_stepping, ECFMetaData()
			.setdescription({ "Auto time stepping" }));

	REGISTER(damping, ECFMetaData()
			.setdescription({ "Damping configuration." }));
}
