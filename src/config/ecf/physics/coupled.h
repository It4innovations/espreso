
#ifndef SRC_CONFIG_ECF_PHYSICS_COUPLED_H_
#define SRC_CONFIG_ECF_PHYSICS_COUPLED_H_

#include "heattransfer.h"
#include "structuralmechanics.h"

namespace espreso {

struct CoupledPhysicsGlobalConfiguration {

	enum class COUPLING {
		WEAK = 0,
		STRONG  = 1
	};

	COUPLING coupling;
};

struct ThermoElasticityLoadStepConfiguration: public CoupledPhysicsGlobalConfiguration, public ECFDescription {

	HeatTransferLoadStepConfiguration heat_transfer;
	StructuralMechanicsLoadStepConfiguration structural_mechanics;

	ThermoElasticityLoadStepConfiguration(DIMENSION *D);
};

struct ThermoElasticityConfiguration: public PhysicsConfiguration, public HeatTransferGlobalSettings, public StructuralMechanicsGlobalSettings {

	std::map<size_t, ThermoElasticityLoadStepConfiguration> load_steps_settings;

	ThermoElasticityConfiguration(DIMENSION dimension);
};

}

#endif /* SRC_CONFIG_ECF_PHYSICS_COUPLED_H_ */
