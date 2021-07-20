
#ifndef SRC_CONFIG_ECF_PHYSICS_ACOUSTIC_H_
#define SRC_CONFIG_ECF_PHYSICS_ACOUSTIC_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct ECF;

struct AcousticGlobalSettings {

	bool init_temp_respect_bc, diffusion_split;

	AcousticGlobalSettings(ECFObject *ecfdescription);
};

struct AcousticLoadStepConfiguration: public AcousticLoadStepSolverConfiguration {

	std::map<std::string, ECFExpression> acoustic_pressure, normal_acceleration;

	AcousticLoadStepConfiguration(DIMENSION *D);
};

struct AcousticConfiguration: public PhysicsConfiguration, public AcousticGlobalSettings {

	DIMENSION dimension;
	std::map<size_t, AcousticLoadStepConfiguration> load_steps_settings;

	AcousticConfiguration(DIMENSION d);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_ACOUSTIC_H_ */
