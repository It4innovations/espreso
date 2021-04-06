
#ifndef SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_

#include "physics.h"
#include "physicssolver/loadstep.h"

namespace espreso {

struct ECF;

struct RotorDynamicsConfiguration: public ECFDescription {
	enum class TYPE {
		FIXED,
		COROTATING
	};

	struct RotationConfiguration: public ECFDescription {
		enum class TYPE {
			FREQUENCY_RATIO,
			TABLE
		};

		TYPE type;
		double frequency_ratio;
		std::vector<double> table;

		RotationConfiguration();
	};

	struct RotationAxisConfiguration: public ECFDescription {
		DIMENSION dimension;
		ECFExpressionVector center, orientation;

		RotationAxisConfiguration(DIMENSION dimension);
	};

	struct CorotatingRotorConfiguration: public ECFDescription {
		bool coriolis_effect, spin_softening, rotating_damping, centrifugal_load;
		RotationConfiguration rotation;

		CorotatingRotorConfiguration();
	};

	struct CorotatingConfiguration: public ECFDescription {
		std::map<std::string, CorotatingRotorConfiguration> rotors_definitions;
		RotationAxisConfiguration rotation_axis;

		CorotatingConfiguration(DIMENSION dimension);
	};

	struct FixedRotorConfiguration: public ECFDescription {
		bool gyroscopic_effect, centrifugal_load;
//		bool rotating_damping;
		RotationConfiguration rotation;
		RotationAxisConfiguration rotation_axis;

		FixedRotorConfiguration(DIMENSION dimension);
	};

	struct FixedConfiguration: public ECFDescription {
		std::map<std::string, FixedRotorConfiguration> rotors_definitions;

		FixedConfiguration(DIMENSION dimension);
	};

	TYPE type;
	FixedConfiguration fixed;
	CorotatingConfiguration corotating;

	RotorDynamicsConfiguration(DIMENSION dimension);
};

struct RotatingForceConfiguration: public ECFDescription {

	ECFExpressionVector rotation_axis;
	ECFExpression rotation_radius;
	ECFExpression unbalance_mass;
	ECFExpression unbalance_phase_angle;
	ECFExpression location;

	RotatingForceConfiguration(DIMENSION *dimension);
};

struct NonlinerSpringConfiguration: public ECFDescription {
	enum class Support {
		FIXED,
		SLIDING
	};

	Support support;
	ECFExpressionVector direction;
	ECFExpression force;
	ECFExpression force_derivative;

	NonlinerSpringConfiguration(DIMENSION *dimension);
};

struct StructuralMechanicsGlobalSettings {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	DIMENSION element_dimension;
	ELEMENT_BEHAVIOUR element_behaviour;

	StructuralMechanicsGlobalSettings(ECFObject *ecfdescription, DIMENSION dimension);
};

struct StructuralMechanicsLoadStepConfiguration: public StructuralMechanicsLoadStepSolverConfiguration {

	bool large_displacement;

	std::map<std::string, ECFExpression> temperature, normal_pressure;
	std::map<std::string, ECFExpressionVector> force, angular_velocity, acceleration, normal_direction, obstacle;
	std::map<std::string, ECFHarmonicExpressionVector> harmonic_force, harmonic_acceleration;
	std::map<std::string, ECFExpressionOptionalVector> displacement;
	std::map<std::string, RotatingForceConfiguration> rotating_force;
	std::map<std::string, NonlinerSpringConfiguration> nonlinear_spring;
	RotorDynamicsConfiguration rotor_dynamics;

	StructuralMechanicsLoadStepConfiguration(DIMENSION *dimension);
};

struct StructuralMechanicsOutputSettings: public virtual ECFDescription {

	static void addMonitorableProperties(ECFMetaData &metadata, const ECF *root);

	bool displacement, stress;

	void basic() {
		displacement = true;
		stress = false;
	}
	void all() {
		displacement = true;
		stress = true;
	}

	static void activate();

	StructuralMechanicsOutputSettings();

protected:
	static bool _activated;
};

struct StructuralMechanicsConfiguration: public PhysicsConfiguration, public StructuralMechanicsGlobalSettings {

	DIMENSION dimension;
	std::map<size_t, StructuralMechanicsLoadStepConfiguration> load_steps_settings;

	StructuralMechanicsConfiguration(DIMENSION d);
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_ */
