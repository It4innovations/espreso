
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

	PhysicsSolver()
	{
		// there is always at least load step 1
		load_steps_settings[1] = new LoadStepSettings<TConvergence>();
	}

	~PhysicsSolver()
	{
		// there is always at least load step 1
		delete load_steps_settings[1];
	}

	SUBMAPTOCONFIG(size_t, LoadStepSettings<TConvergence>, load_steps_settings, "Detail settings for each load step.", "1", "Description of load step '1'");
};

}




#endif /* SRC_CONFIGURATION_PHYSICS_SOLVER_H_ */
