
#ifndef SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_DAMPING_H_
#define SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_DAMPING_H_

#include "transientfirstorderimplicit.h"
#include "config/holders/expression.h"

namespace espreso {

struct DirectDampingConfiguration: public ECFDescription {

	ECFExpression mass, stiffness;

	DirectDampingConfiguration();
};

struct RatioDampingConfiguration: public ECFDescription {

	ECFExpression ratio, frequency;

	RatioDampingConfiguration();
};

struct RayleighDampingConfiguration: public ECFDescription {

	enum class TYPE {
		NONE,
		DIRECT,
		DAMPING_RATIO
	};

	TYPE type;

	DirectDampingConfiguration direct_damping;
	RatioDampingConfiguration ratio_damping;

	RayleighDampingConfiguration();
};

struct HarmonicRayleighDampingConfiguration: public RayleighDampingConfiguration {

	ECFExpression structural_damping_coefficient;

	HarmonicRayleighDampingConfiguration();
};

struct CoriolisEffectAxisConfiguration: public ECFDescription {

	double x, y, z;

	CoriolisEffectAxisConfiguration();
};

struct CoriolisEffectConfiguration: public ECFDescription {

	bool coriolis_damping;
	bool spin_softening;

	CoriolisEffectAxisConfiguration rotation_axis;

	CoriolisEffectConfiguration();
};

struct DampingConfiguration: public ECFDescription {

	RayleighDampingConfiguration rayleigh;
//	CoriolisEffectConfiguration coriolis_effect;

	DampingConfiguration();
};

struct HarmonicDampingConfiguration: public ECFDescription {

	HarmonicRayleighDampingConfiguration rayleigh;
//	CoriolisEffectConfiguration coriolis_effect;

	HarmonicDampingConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_PHYSICS_PHYSICSSOLVER_DAMPING_H_ */
