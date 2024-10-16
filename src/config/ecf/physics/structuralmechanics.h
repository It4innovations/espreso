
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
		ECFExpressionVector center, orientation;

		RotationAxisConfiguration();
	};

	struct CorotatingRotorConfiguration: public ECFDescription {
		bool coriolis_effect, spin_softening, rotating_damping, centrifugal_load;
		RotationConfiguration rotation;

		CorotatingRotorConfiguration();
	};

	struct CorotatingConfiguration: public ECFDescription {
		std::map<std::string, CorotatingRotorConfiguration> rotors_definitions;
		RotationAxisConfiguration rotation_axis;

		CorotatingConfiguration();
	};

	struct FixedRotorConfiguration: public ECFDescription {
		bool gyroscopic_effect, centrifugal_load;
//		bool rotating_damping;
		RotationConfiguration rotation;
		RotationAxisConfiguration rotation_axis;

		FixedRotorConfiguration();
	};

	struct FixedConfiguration: public ECFDescription {
		std::map<std::string, FixedRotorConfiguration> rotors_definitions;

		FixedConfiguration();
	};

	TYPE type;
	FixedConfiguration fixed;
	CorotatingConfiguration corotating;

	RotorDynamicsConfiguration();
};

struct RotatingForceConfiguration: public ECFDescription {

	ECFExpressionVector rotation_axis;
	ECFExpression rotation_radius;
	ECFExpression unbalance_mass;
	ECFExpression unbalance_phase_angle;
	ECFExpression location;

	RotatingForceConfiguration();
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

	NonlinerSpringConfiguration();
};

struct FixedWallConfiguration: public ECFDescription {
	ECFExpressionVector normal, point;
	double gap;

	FixedWallConfiguration();
};

struct PressureConfiguration: public ECFDescription {
    ECFExpression pressure;
    ECFExpressionVector direction;

    PressureConfiguration();
};

struct StructuralMechanicsGlobalSettings {

	enum class ELEMENT_BEHAVIOUR {
		PLANE_STRAIN = 0,
		AXISYMMETRIC = 1,
		PLANE_STRESS = 2,
		PLANE_STRESS_WITH_THICKNESS = 3
	};

	ELEMENT_BEHAVIOUR element_behaviour;

	std::map<std::string, ECFExpressionVector> initial_velocity;

	StructuralMechanicsGlobalSettings(ECFObject *ecfdescription);
};

struct StructuralMechanicsLoadStepConfiguration: public StructuralMechanicsLoadStepSolverConfiguration {

	bool large_displacement, corotation;

	std::map<std::string, ECFExpression> temperature, normal_pressure;
	std::map<std::string, PressureConfiguration> pressure;
	std::map<std::string, ECFExpressionVector> force, angular_velocity, acceleration;
	std::map<std::string, FixedWallConfiguration> fixed_wall;
	std::map<std::string, ECFHarmonicExpressionVector> harmonic_force, harmonic_acceleration;
	std::map<std::string, ECFExpressionOptionalVector> displacement;
	std::map<std::string, RotatingForceConfiguration> rotating_force;
	std::map<std::string, NonlinerSpringConfiguration> nonlinear_spring;
	RotorDynamicsConfiguration rotor_dynamics;

	StructuralMechanicsLoadStepConfiguration();
};

struct StructuralMechanicsOutputSettings: public virtual ECFDescription {

	static void addMonitorableProperties(ECFMetaData &metadata, const ECF *root);

	bool displacement, stress, reactions;

	void basic() {
		displacement = true;
		stress = false;
		reactions = false;
	}
	void all() {
		displacement = true;
		stress = true;
		reactions = true;
	}

	static void activate();

	StructuralMechanicsOutputSettings();

protected:
	static bool _activated;
};

struct StructuralMechanicsConfiguration: public PhysicsConfiguration, public StructuralMechanicsGlobalSettings {

	std::map<size_t, StructuralMechanicsLoadStepConfiguration> load_steps_settings;

	StructuralMechanicsConfiguration();
};

}


#endif /* SRC_CONFIG_ECF_PHYSICS_STRUCTURALMECHANICS_H_ */
