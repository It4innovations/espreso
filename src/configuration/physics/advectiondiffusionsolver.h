
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_

#include "nonlinearsolver.h"
#include "transientsolver.h"

namespace espreso {

struct AdvectionDiffusionNonLinearConvergence: public NonLinearConvergence {

	virtual bool checkSolution() const { return temperature; }
	virtual bool checkResidual() const { return heat; }

	virtual double requestedSolution() const { return temperature_residual; }
	virtual double requestedResidual() const { return heat_residual; }

	PARAMETER(bool, temperature, "Turn on/off temperature residual check.", true);
	PARAMETER(bool, heat       , "Turn on/off heat residual check."       , false);

	PARAMETER(double, temperature_residual, "Requested temperature residual", 1e-3);
	PARAMETER(double, heat_residual       , "Requested heat residual"       , 1e-3);
};

struct AdvectionDiffusionSolver: public Configuration {

	enum class TYPE {
		STEADY_STATE,
		TRANSIENT
	};

	enum class MODE {
		LINEAR,
		NONLINEAR
	};

	enum class INTERPOLATION {
		LINEAR,
		QUADRATIC
	};

	OPTION(TYPE, type, "Solver type", TYPE::STEADY_STATE, OPTIONS({
		{ "STEADY_STATE", TYPE::STEADY_STATE, "Steady state." },
		{ "TRANSIENT"   , TYPE::TRANSIENT   , "Transient." },
	}));

	OPTION(MODE, mode, "Solver mode", MODE::LINEAR, OPTIONS({
		{ "LINEAR"   , MODE::LINEAR   , "Linear." },
		{ "NONLINEAR", MODE::NONLINEAR, "Non-linear." },
	}));

	OPTION(INTERPOLATION, interpolation, "Interpolation of tabular data.", INTERPOLATION::LINEAR, OPTIONS({
		{ "LINEAR"   , INTERPOLATION::LINEAR   , "Linear." },
		{ "QUADRATIC", INTERPOLATION::QUADRATIC, "Quadratic." },
	}));

	PARAMETER(size_t, load_steps, "Number of load steps in simulation.", 1);

	SUBMAPTOCONFIG(size_t, NonLinearSolver<AdvectionDiffusionNonLinearConvergence>, nonlinear_solver, "Non-linear configuration for each load step.", "1", "Description of non-linear solver for step '1'");
	SUBMAPTOCONFIG(size_t, TransientSolver, transient_solver, "Transient configuration for each load step.", "1", "Description of transient solver for step '1'");
};

}


#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_ */
