
#ifndef SRC_CONFIG_ECF_PHYSICS_ACOUSTICS_H_
#define SRC_CONFIG_ECF_PHYSICS_ACOUSTICS_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct ECF;

struct AcousticsLoadStepConfiguration: public LoadStepSolverConfiguration {

	std::map<std::string, ECFExpression> acoustic_pressure, normal_acceleration;

	AcousticsLoadStepConfiguration(DIMENSION *D);
};

struct AcousticsConfiguration: public PhysicsConfiguration {

	DIMENSION dimension;
	std::map<size_t, AcousticsLoadStepConfiguration> load_steps_settings;

	AcousticsConfiguration(DIMENSION d);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_ACOUSTICS_H_ */
