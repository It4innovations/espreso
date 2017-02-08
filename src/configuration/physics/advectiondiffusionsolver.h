
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_

#include "nonlinearsolver.h"
#include "transientsolver.h"

namespace espreso {

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

	SUBMAPTOCONFIG(size_t, NonLinearSolver, nonlinear_solver, "Non-linear configuration for each load step.");
	SUBMAPTOCONFIG(size_t, TransientSolver, transient_solver, "Transient configuration for each load step.");
};

}


#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONSOLVER_H_ */
