
#ifndef SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_

#include "../configuration.h"

namespace espreso {

struct NonLinearConvergence: public Configuration {

	PARAMETER(bool, temperature, "Turn on/off temperature residual check.", true);
	PARAMETER(bool, heat       , "Turn on/off heat residual check."       , false);

	PARAMETER(double, temperature_residual, "Requested temperature residual", 1e-3);
	PARAMETER(double, heat_residual       , "Requested heat residual"       , 1e-3);
};

struct NonLinearSolver: public Configuration {

	enum class METHOD {
		NEWTON_RHAPSON
	};

	OPTION(METHOD, method, "Non-linear method", METHOD::NEWTON_RHAPSON, OPTIONS({
		{ "NEWTON_RHAPSON", METHOD::NEWTON_RHAPSON, "Newton-Rhapson method." }
	}));

	PARAMETER(size_t, max_iterations, "Allowed number of iterations.", 100);
	PARAMETER(bool, line_search, "Set line search.", false);
	SUBCONFIG(NonLinearConvergence, convergence_parameters, "Convergence parameters.");
	PARAMETER(size_t, substeps, "Number of loading substeps.", 1);
};

}

#endif /* SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_ */
