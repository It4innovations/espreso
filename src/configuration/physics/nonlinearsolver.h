
#ifndef SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_

#include "../configuration.hpp"

namespace espreso {

struct NonLinearConvergence: public Configuration {

	virtual bool checkSolution() const =0;
	virtual bool checkResidual() const =0;

	virtual double requestedSolution() const =0;
	virtual double requestedResidual() const =0;
};

struct NonLinearSolverBase: public Configuration {

	virtual const NonLinearConvergence& convergenceParameters() const =0;

	enum class METHOD {
		NEWTON_RHAPSON,
		MODIFIED_NEWTON_RHAPSON
	};

	enum class STEPPINGG {
		TRUE,
		FALSE,
		AUTO
	};

	OPTION(METHOD, method, "Non-linear method", METHOD::NEWTON_RHAPSON, OPTIONS({
		{ "NEWTON_RHAPSON"         , METHOD::NEWTON_RHAPSON         , "Newton-Rhapson method." },
		{ "MODIFIED_NEWTON_RHAPSON", METHOD::MODIFIED_NEWTON_RHAPSON, "Modified Newton-Rhapson method (constant stifness matrix)." }
	}));

	OPTION(STEPPINGG, stepping, "Non-linear sub-stepping", STEPPINGG::FALSE, OPTIONS({
		{ "FALSE", STEPPINGG::FALSE, "Turn off stepping." },
		{ "TRUE" , STEPPINGG::TRUE , "Turn on stepping." },
		{ "AUTO" , STEPPINGG::AUTO , "Automatic stepping." }
	}));

	PARAMETER(size_t, max_iterations, "Allowed number of iterations.", 15);
	PARAMETER(bool, line_search, "Set line search.", false);
	PARAMETER(size_t, substeps, "Number of loading substeps.", 1);
};


template <class TConvergence>
struct NonLinearSolver: public NonLinearSolverBase {

	const NonLinearConvergence& convergenceParameters() const { return convergence_parameters; }

	SUBCONFIG(TConvergence, convergence_parameters, "Convergence parameters.");
};

}

#endif /* SRC_CONFIGURATION_PHYSICS_NONLINEARSOLVER_H_ */
