
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

struct ImpedanceConfiguration: public ECFDescription {
	ECFExpression impedance;
	ImpedanceConfiguration();
};

struct PointSourceConfiguration: public ECFDescription {

	enum class TYPE {
		USER,
		POWER,
		INTENSITY,
		FLOW
	};

	TYPE type;
	ECFExpression volume_flow, phase, reference_intensity, distance_from_source, reference_power, source_amplitude;

	PointSourceConfiguration();
};

struct AcousticLoadStepConfiguration: public AcousticLoadStepSolverConfiguration {

	enum class SYSTEM {
		REAL,
		COMPLEX
	};

	SYSTEM system;

	std::map<std::string, ECFExpression> acoustic_pressure, normal_acceleration, monopole_source;
	std::map<std::string, ECFExpressionVector> dipole_source, acceleration;
	std::map<std::string, ImpedanceConfiguration> impedance;
	std::map<std::string, PointSourceConfiguration> point_source;

	AcousticLoadStepConfiguration(DIMENSION *D);
};

struct AcousticConfiguration: public PhysicsConfiguration, public AcousticGlobalSettings {

	DIMENSION dimension;
	std::map<size_t, AcousticLoadStepConfiguration> load_steps_settings;

	AcousticConfiguration(DIMENSION d);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_ACOUSTIC_H_ */
