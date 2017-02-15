
#ifndef SRC_CONFIGURATION_PHYSICS_SOLVER_H_
#define SRC_CONFIGURATION_PHYSICS_SOLVER_H_

#include "loadstepsettings.h"

namespace espreso {

struct PhysicsSolverBase: public Configuration {

	enum class INTERPOLATION {
		LINEAR,
		QUADRATIC
	};

	OPTION(INTERPOLATION, interpolation, "Interpolation of tabular data.", INTERPOLATION::LINEAR, OPTIONS({
		{ "LINEAR"   , INTERPOLATION::LINEAR   , "Linear." },
		{ "QUADRATIC", INTERPOLATION::QUADRATIC, "Quadratic." },
	}));

	PARAMETER(size_t, load_steps, "Number of load steps in simulation.", 1);
};

template <class TConvergence>
struct PhysicsSolver: public PhysicsSolverBase {

	SUBMAPTOCONFIG(size_t, LoadStepSettings<TConvergence>, load_steps_settings, "Detail settings for each load step.");
};

}




#endif /* SRC_CONFIGURATION_PHYSICS_SOLVER_H_ */
