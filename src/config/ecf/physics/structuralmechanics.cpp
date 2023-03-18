
#include "config/ecf/ecf.h"
#include "config/configuration.hpp"

using namespace espreso;

RotorDynamicsConfiguration::RotationAxisConfiguration::RotationAxisConfiguration(DIMENSION dimension)
: dimension(dimension), center(&this->dimension), orientation(&this->dimension)
{
	REGISTER(center, ECFMetaData()
			.setname("Center")
			.setdescription({ "A rotation center of the rotor." }));

	REGISTER(orientation, ECFMetaData()
			.setname("Orientation.")
			.setdescription({ "An orientation of the rotor." }));
}

RotorDynamicsConfiguration::RotationConfiguration::RotationConfiguration()
{
	type = TYPE::FREQUENCY_RATIO;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Rotation configuration." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FREQUENCY_RATIO").setdescription("Frequency ration."))
			.addoption(ECFOption().setname("TABLE").setdescription("A list of frequencies.")));

	frequency_ratio = 1;
	REGISTER(frequency_ratio, ECFMetaData()
			.setdescription({ "Frequency ratio" })
			.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(table, ECFMetaData()
			.setdescription({ "Requested frequencies." })
			.setdatatype({ ECFDataType::STRING }));
}

RotorDynamicsConfiguration::CorotatingRotorConfiguration::CorotatingRotorConfiguration()
{
	coriolis_effect = true;
	REGISTER(coriolis_effect, ECFMetaData()
			.setdescription({ "Coriolis effect" })
			.setdatatype({ ECFDataType::BOOL }));

	spin_softening = true;
	REGISTER(spin_softening, ECFMetaData()
			.setdescription({ "Spin softening" })
			.setdatatype({ ECFDataType::BOOL }));

	rotating_damping = false;
//	REGISTER(rotating_damping, ECFMetaData()
//			.setdescription({ "Rotation damping" })
//			.setdatatype({ ECFDataType::BOOL }));

	centrifugal_load = true;
	REGISTER(centrifugal_load, ECFMetaData()
			.setdescription({ "centrifugal load" })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(rotation, ECFMetaData().setdescription({ "Rotation configuration." }));
}

RotorDynamicsConfiguration::CorotatingConfiguration::CorotatingConfiguration(DIMENSION dimension)
: rotation_axis(dimension)
{
	REGISTER(rotors_definitions, ECFMetaData()
			.setdescription({ "List of rotor definitions.", "Rotor definition." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }));

	REGISTER(rotation_axis, ECFMetaData().setdescription({ "Rotation axis." }));
}

RotorDynamicsConfiguration::FixedRotorConfiguration::FixedRotorConfiguration(DIMENSION dimension)
: rotation_axis(dimension)
{
	gyroscopic_effect = true;
	REGISTER(gyroscopic_effect, ECFMetaData()
			.setdescription({ "Rotation damping" })
			.setdatatype({ ECFDataType::BOOL }));

//	rotating_damping = true;
//	REGISTER(rotating_damping, ECFMetaData()
//			.setdescription({ "Rotation damping" })
//			.setdatatype({ ECFDataType::BOOL }));

	centrifugal_load = true;
	REGISTER(centrifugal_load, ECFMetaData()
			.setdescription({ "Rotation damping" })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(rotation, ECFMetaData().setdescription({ "Rotation configuration." }));
	REGISTER(rotation_axis, ECFMetaData().setdescription({ "Rotation axis." }));
}

RotorDynamicsConfiguration::FixedConfiguration::FixedConfiguration(DIMENSION dimension)
{
	REGISTER(rotors_definitions, ECFMetaData()
			.setdescription({ "List of rotor definitions.", "Rotor definition." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY REGIONS" }), dimension);
}

RotorDynamicsConfiguration::RotorDynamicsConfiguration(DIMENSION dimension)
: fixed(dimension), corotating(dimension)
{
	type = TYPE::FIXED;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Rotor type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FIXED").setdescription("Fixed rotor."))
			.addoption(ECFOption().setname("COROTATING").setdescription("Corotating rotor.")));

	REGISTER(fixed, ECFMetaData().setdescription({ "Fixed rotor configuration." }));
	REGISTER(corotating, ECFMetaData().setdescription({ "Corotating rotor configuration." }));
}

RotatingForceConfiguration::RotatingForceConfiguration(DIMENSION *dimension)
: rotation_axis(dimension)
{
	rotation_axis.x.value = rotation_axis.y.value = "0";
	rotation_axis.z.value = "1";
	REGISTER(rotation_axis, ECFMetaData().setdescription({ "Rotation axis." }));
	rotation_radius.value = "0";
	REGISTER(rotation_radius, ECFMetaData()
			.setdescription({ "Rotation radius." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	unbalance_mass.value = "0";
	REGISTER(unbalance_mass, ECFMetaData()
			.setdescription({ "Unbalance mass." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	unbalance_phase_angle.value = "0";
	REGISTER(unbalance_phase_angle, ECFMetaData()
			.setdescription({ "Unbalance phase angle." })
			.setdatatype({ ECFDataType::EXPRESSION }));
	location.value = "0";
	REGISTER(location, ECFMetaData()
			.setdescription({ "location." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

NonlinerSpringConfiguration::NonlinerSpringConfiguration(DIMENSION *dimension)
: direction(dimension)
{
	support = Support::SLIDING;
	REGISTER(support, ECFMetaData()
			.setdescription({ "Support type." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("FIXED").setdescription("Fixed."))
			.addoption(ECFOption().setname("SLIDING").setdescription("Sliding.")));

	REGISTER(direction, ECFMetaData().setdescription({ "Direction." }));
	REGISTER(force, ECFMetaData().setdescription({ "Force." }).setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(force_derivative, ECFMetaData().setdescription({ "Force derivative." }).setdatatype({ ECFDataType::EXPRESSION }));
}

StructuralMechanicsLoadStepConfiguration::StructuralMechanicsLoadStepConfiguration(DIMENSION *dimension)
: rotor_dynamics(*dimension)
{
	large_displacement = false;
	REGISTER(large_displacement, ECFMetaData()
			.setdescription({ "Turn on large displacement." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Temperature of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "275.15" }));
	REGISTER(normal_pressure, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal pressure on a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "0" }));

	REGISTER(force, ECFMetaData()
			.setdescription({ "The name of a region.", "Force on a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);

	REGISTER(angular_velocity, ECFMetaData()
			.setdescription({ "The name of a region.", "Angular velocity of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);

	REGISTER(normal_direction, ECFMetaData()
			.setdescription({ "The name of a region.", "Normal direction of a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);

	REGISTER(obstacle, ECFMetaData()
			.setdescription({ "The name of a region.", "Obstacle for a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);

	REGISTER(acceleration, ECFMetaData()
			.setdescription({ "The name of a region.", "Acceleration of a given region." })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);


	REGISTER(displacement, ECFMetaData()
			.setdescription({ "The name of a region.", "Fixed displacement of a given region." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" }),
			dimension);

	REGISTER(harmonic_force, ECFMetaData()
			.setdescription({ "The name of a region.", "Harmonic force" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION", }),
			dimension);

	REGISTER(nonlinear_spring, ECFMetaData()
			.setdescription({ "The name of a region.", "Non-linear spring" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION", }),
			dimension);

	REGISTER(rotating_force, ECFMetaData()
			.setdescription({ "The name of a region.", "Rotating force" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION", }),
			dimension);

	REGISTER(rotor_dynamics, ECFMetaData().setdescription({ "Rotor dynamics." }));
}

bool StructuralMechanicsOutputSettings::_activated = false;

void StructuralMechanicsOutputSettings::addMonitorableProperties(ECFMetaData &metadata, const ECF *root)
{
	auto harmonic = [&] () {
		if (root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
			for (auto it = root->structural_mechanics_3d.load_steps_settings.begin(); it != root->structural_mechanics_3d.load_steps_settings.end(); ++it) {
				if (it->second.type == LoadStepSolverConfiguration::TYPE::HARMONIC) {
					return true;
				}
			}
		}
		return false;
	};

	metadata
	.addoption(ECFOption().setname("DISPLACEMENT").setdescription("Displacement magnitude.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("DISPLACEMENT_X").setdescription("Displacement in x-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("DISPLACEMENT_Y").setdescription("Displacement in y-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("DISPLACEMENT_Z").setdescription("Displacement in z-direction.").allowonly([&] () { return _activated && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("DISPLACEMENT_AMPLITUDE_X").setdescription("Displacement amplitude in x-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_AMPLITUDE_Y").setdescription("Displacement amplitude in y-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_AMPLITUDE_Z").setdescription("Displacement amplitude in z-direction.").allowonly([&] () { return _activated && harmonic() && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("DISPLACEMENT_COS_X").setdescription("Displacement cosine in x-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_COS_Y").setdescription("Displacement cosine in y-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_COS_Z").setdescription("Displacement cosine in z-direction.").allowonly([&] () { return _activated && harmonic() && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("DISPLACEMENT_SIN_X").setdescription("Displacement sinus in x-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_SIN_Y").setdescription("Displacement sinus in y-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("DISPLACEMENT_SIN_Z").setdescription("Displacement sinus in z-direction.").allowonly([&] () { return _activated && harmonic() && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("VELOCITY").setdescription("Velocity magnitude.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("VELOCITY_X").setdescription("Velocity in x-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("VELOCITY_Y").setdescription("Velocity in y-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("VELOCITY_Z").setdescription("Velocity in z-direction.").allowonly([&] () { return _activated && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("VELOCITY_AMPLITUDE_X").setdescription("Velocity amplitude in x-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("VELOCITY_AMPLITUDE_Y").setdescription("Velocity amplitude in y-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("VELOCITY_AMPLITUDE_Z").setdescription("Velocity amplitude in z-direction.").allowonly([&] () { return _activated && harmonic() && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("ACCELERATION").setdescription("Acceleration magnitude.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("ACCELERATION_X").setdescription("Acceleration in x-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("ACCELERATION_Y").setdescription("Acceleration in y-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("ACCELERATION_Z").setdescription("Acceleration in z-direction.").allowonly([&] () { return _activated && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }))

	.addoption(ECFOption().setname("ACCELERATION_AMPLITUDE_X").setdescription("Acceleration amplitude in x-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("ACCELERATION_AMPLITUDE_Y").setdescription("Acceleration amplitude in y-direction.").allowonly([&] () { return _activated && harmonic(); }))
	.addoption(ECFOption().setname("ACCELERATION_AMPLITUDE_Z").setdescription("Acceleration amplitude in z-direction.").allowonly([&] () { return _activated && harmonic() && root->physics == PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D; }));
}


StructuralMechanicsOutputSettings::StructuralMechanicsOutputSettings()
{
	basic();

	REGISTER(displacement, ECFMetaData()
				.setdescription({ "Displacement." })
				.setdatatype({ ECFDataType::BOOL })
				.allowonly([&] () { return _activated; }));

	REGISTER(stress, ECFMetaData()
				.setdescription({ "Stress." })
				.setdatatype({ ECFDataType::BOOL })
				.allowonly([&] () { return _activated; }));
}

StructuralMechanicsGlobalSettings::StructuralMechanicsGlobalSettings(ECFObject *ecfdescription, DIMENSION dimension)
: element_dimension(dimension)
{
	element_behaviour = ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS;
	REGISTER(element_behaviour, ECFMetaData()
			.setdescription({ "Physics solver type." })
			.setdatatype({ ECFDataType::OPTION })
			.allowonly([&] () { return element_dimension == DIMENSION::D2; })
			.addoption(ECFOption().setname("PLANE_STRAIN").setdescription("Plane strain."))
			.addoption(ECFOption().setname("AXISYMMETRIC").setdescription("Axisymmetric."))
			.addoption(ECFOption().setname("PLANE_STRESS").setdescription("Plane stress."))
			.addoption(ECFOption().setname("PLANE_STRESS_WITH_THICKNESS").setdescription("Plane stress with thickness.")));
}

StructuralMechanicsConfiguration::StructuralMechanicsConfiguration(DIMENSION d)
: PhysicsConfiguration(d, MaterialConfiguration::PHYSICAL_MODEL::STRUCTURAL_MECHANICS),
  StructuralMechanicsGlobalSettings(ecfdescription, d),
  dimension(d)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step.", "Load step index." })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" }),
			&dimension);
}




