
#ifndef SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_
#define SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_

#include "nonlinearsolver.h"
#include "transientsolver.h"

namespace espreso {

struct LoadStepSettings: public Configuration {

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

	PARAMETER(size_t, start_time, "Step start time.", 0);
	PARAMETER(size_t, end_time  , "Step end time"   , 1);

	SUBCONFIG(NonLinearSolver, nonlinear_solver, "Non-linear configuration for each load step.");
	SUBCONFIG(TransientSolver, transient_solver, "Transient configuration for each load step.");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_LOADSTEPSETTINGS_H_ */
