
#ifndef SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_
#define SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_

#include "nonlinearsolver.h"
#include "transientsolver.h"
#include "../solver.h"

namespace espreso {

struct LoadStepSettingsBase: public Configuration {

	enum class TYPE {
		STEADY_STATE,
		TRANSIENT
	};

	enum class MODE {
		LINEAR,
		NONLINEAR
	};

	OPTION(TYPE, type, "Solver type", TYPE::STEADY_STATE, OPTIONS({
		{ "STEADY_STATE", TYPE::STEADY_STATE, "Steady state." },
		{ "TRANSIENT"   , TYPE::TRANSIENT   , "Transient." },
	}));

	OPTION(MODE, mode, "Solver mode", MODE::LINEAR, OPTIONS({
		{ "LINEAR"   , MODE::LINEAR   , "Linear." },
		{ "NONLINEAR", MODE::NONLINEAR, "Non-linear." },
	}));

	OPTION(SOLVER_LIBRARY, solver_library, "Linear solver used for computing a system.", SOLVER_LIBRARY::ESPRESO, OPTIONS({
		{ "ESPRESO", SOLVER_LIBRARY::ESPRESO, "ESPRESO solver [FETI methods]" },
		{ "HYPRE"  , SOLVER_LIBRARY::HYPRE  , "Hypre solver [multigrid methods]" },
	}));

	SUBCONFIG(ESPRESOSolver, espreso, "Internal FETI solver options.");
	SUBCONFIG(HypreSolver  , hypre  , "Multigrid solver setting.");

	PARAMETER(double, duration_time, "Duration time of the load step.", 1);

	SUBCONFIG(TransientSolver, transient_solver, "Transient configuration for each load step.");
};

template<class TConvergence>
struct LoadStepSettings: public LoadStepSettingsBase {

	SUBCONFIG(NonLinearSolver<TConvergence>, nonlinear_solver, "Transient configuration for each load step.");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_ */
