
#include <config/ecf/ecf.h>
#include "config/configuration.hpp"

using namespace espreso;

ConvectionConfiguration::ConvectionConfiguration()
: heat_transfer_coefficient(ECFMetaData::getboundaryconditionvariables()),
  external_temperature(ECFMetaData::getboundaryconditionvariables()),
  wall_height(ECFMetaData::getboundaryconditionvariables()),
  tilt_angle(ECFMetaData::getboundaryconditionvariables()),
  diameter(ECFMetaData::getboundaryconditionvariables()),
  plate_length(ECFMetaData::getboundaryconditionvariables()),
  fluid_velocity(ECFMetaData::getboundaryconditionvariables()),
  plate_distance(ECFMetaData::getboundaryconditionvariables()),
  length(ECFMetaData::getboundaryconditionvariables()),
  experimental_constant(ECFMetaData::getboundaryconditionvariables()),
  volume_fraction(ECFMetaData::getboundaryconditionvariables()),
  absolute_pressure(ECFMetaData::getboundaryconditionvariables())
{
	type = TYPE::USER;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Convection type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("USER").setdescription("User defined."))
			.addoption(ECFOption().setname("EXTERNAL_NATURAL").setdescription("External natural."))
			.addoption(ECFOption().setname("INTERNAL_NATURAL").setdescription("Internal natural."))
			.addoption(ECFOption().setname("EXTERNAL_FORCED").setdescription("External forced."))
			.addoption(ECFOption().setname("INTERNAL_FORCED").setdescription("Internal forced.")));

	variant = VARIANT::VERTICAL_WALL;
	REGISTER(variant, ECFMetaData()
			.setdescription({ "Variant" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("VERTICAL_WALL").setdescription("Vertical wall."))
			.addoption(ECFOption().setname("INCLINED_WALL").setdescription("Inclined wall."))
			.addoption(ECFOption().setname("HORIZONTAL_CYLINDER").setdescription("Horizontal cylinder."))
			.addoption(ECFOption().setname("SPHERE").setdescription("Sphere."))
			.addoption(ECFOption().setname("HORIZONTAL_PLATE_UP").setdescription("Horizontal place up."))
			.addoption(ECFOption().setname("HORIZONTAL_PLATE_DOWN").setdescription("Horizontal plate down."))
			.addoption(ECFOption().setname("AVERAGE_PLATE").setdescription("Average plate."))
			.addoption(ECFOption().setname("PARALLEL_PLATES").setdescription("Parallel plates."))
			.addoption(ECFOption().setname("CIRCULAR_TUBE").setdescription("Circular tube."))
			.addoption(ECFOption().setname("TUBE").setdescription("Tube."))
			.addoption(ECFOption().setname("QUENCH_PARALLEL").setdescription("Quench parallel plate."))
			.addoption(ECFOption().setname("QUENCH").setdescription("Quench.")));

	fluid = FLUID::AIR;
	REGISTER(fluid, ECFMetaData()
			.setdescription({ "Fluid type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("AIR").setdescription("Air."))
			.addoption(ECFOption().setname("WATER").setdescription("Water."))
			.addoption(ECFOption().setname("ENGINE_OIL").setdescription("Engine oil."))
			.addoption(ECFOption().setname("TRANSFORMER_OIL").setdescription("Tranformer oil."))
			.addoption(ECFOption().setname("STEAM").setdescription("Steam.")));

	REGISTER(heat_transfer_coefficient, ECFMetaData()
			.setdescription({ "Heat transfer coefficient" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(external_temperature, ECFMetaData()
			.setdescription({ "Ambient temperature" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(wall_height, ECFMetaData()
			.setdescription({ "Wall height" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(tilt_angle, ECFMetaData()
			.setdescription({ "Tilt angle" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(diameter, ECFMetaData()
			.setdescription({ "Diameter" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(plate_length, ECFMetaData()
			.setdescription({ "Plate length" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(fluid_velocity, ECFMetaData()
			.setdescription({ "Fluid velocity" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(plate_distance, ECFMetaData()
			.setdescription({ "Plate distance" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(length, ECFMetaData()
			.setdescription({ "Length" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(experimental_constant, ECFMetaData()
			.setdescription({ "Experimental constant" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(volume_fraction, ECFMetaData()
			.setdescription({ "Warer volume fraction" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(absolute_pressure, ECFMetaData()
			.setdescription({ "Absolute pressure" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

RadiationConfiguration::RadiationConfiguration()
: emissivity(ECFMetaData::getboundaryconditionvariables()),
  external_temperature(ECFMetaData::getboundaryconditionvariables())
{
	REGISTER(emissivity, ECFMetaData()
			.setdescription({ "Emissivity" })
			.setdatatype({ ECFDataType::EXPRESSION }));
	REGISTER(external_temperature, ECFMetaData()
			.setdescription({ "Ambient temperature" })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

BioHeatSourceConfiguration::BioHeatSourceConfiguration()
: arteriar_blood_temperature(ECFMetaData::getboundaryconditionvariables()),
  blood_specific_heat(ECFMetaData::getboundaryconditionvariables()),
  blood_density(ECFMetaData::getboundaryconditionvariables()),
  metabolic_heat_source(ECFMetaData::getboundaryconditionvariables()),
  blood_perfusion(ECFMetaData::getboundaryconditionvariables()),
  reference_temperature(ECFMetaData::getboundaryconditionvariables()),
  physical_activity_scatter_factor(ECFMetaData::getboundaryconditionvariables()),
  mu(ECFMetaData::getboundaryconditionvariables())
{
	REGISTER(arteriar_blood_temperature, ECFMetaData()
		.setdescription({ "Arteriar blood temperature." })
		.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(blood_specific_heat, ECFMetaData()
			.setdescription({ "Blood specific heat." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(blood_density, ECFMetaData()
			.setdescription({ "Blood density." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(metabolic_heat_source, ECFMetaData()
			.setdescription({ "Metabolic heat source." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(blood_perfusion, ECFMetaData()
			.setdescription({ "Blood perfusion." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(reference_temperature, ECFMetaData()
			.setdescription({ "Reference temperature." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(physical_activity_scatter_factor, ECFMetaData()
			.setdescription({ "scatter factor of muscle work do to physical activity." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(mu, ECFMetaData()
			.setdescription({ "MU." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

HumanThermoregulationSystem::HumanThermoregulationSystem()
{
	activity_level_unit = ACTIVITY_LEVEL_UNIT::STANDING;
	REGISTER(activity_level_unit, ECFMetaData()
			.setdescription({ "human thermoregulation system settings" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("STANDING").setdescription("--."))
			.addoption(ECFOption().setname("EASY_WORK").setdescription("--."))
			.addoption(ECFOption().setname("HARD_WORK").setdescription("--."))
			.addoption(ECFOption().setname("WALKING_0_PRCT").setdescription("--."))
			.addoption(ECFOption().setname("WALKING_5_PRCT").setdescription("--."))
			.addoption(ECFOption().setname("WALKING_15_PRCT").setdescription("--."))
			.addoption(ECFOption().setname("TEACHER").setdescription("--")));
}

HeatTransferLoadStepConfiguration::HeatTransferLoadStepConfiguration(DIMENSION *D)
{
	update_initial_temperature = false;
	REGISTER(update_initial_temperature, ECFMetaData()
			.setdescription({ "Set initial temperature as the output from previous load step." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "The name of a region.", "Temperature" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(heat_source, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat source" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "273.15" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(translation_motions, ECFMetaData()
			.setdescription({ "The name of a region.", "Translation motion" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" })
			.setdynamic(),
			D, ECFMetaData::getboundaryconditionvariables(), "0");

	REGISTER(heat_flux, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat flux" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "500" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());
	REGISTER(heat_flow, ECFMetaData()
			.setdescription({ "The name of a region.", "Heat flow" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION, ECFDataType::EXPRESSION })
			.setpattern({ "MY_REGION", "500" })
			.setdynamic(),
			ECFMetaData::getboundaryconditionvariables());

	REGISTER(bio_heat, ECFMetaData()
			.setdescription({ "The name of a region.", "Bio heat" })
			.setdatatype({ ECFDataType::ELEMENTS_REGION })
			.setpattern({ "MY_REGION" })
			.setdynamic());

	REGISTER(convection, ECFMetaData()
			.setdescription({ "The name of a region.", "Convection" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" })
			.setdynamic());
	REGISTER(diffuse_radiation, ECFMetaData()
			.setdescription({ "The name of a region.", "Diffuse radiation" })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "MY_REGION" })
			.setdynamic());

	REGISTER(human_thermoregulation_system, ECFMetaData()
			.setdescription({  "human thermoregulation system settings" })
			.setcollapsed());

}

bool HeatTransferOutputSettings::_activated = false;

void HeatTransferOutputSettings::addMonitorableProperties(ECFMetaData &metadata, const ECF *root)
{
	metadata
	.addoption(ECFOption().setname("TEMPERATURE").setdescription("Minimum.").allowonly([&] () { return _activated; }))

	.addoption(ECFOption().setname("FLUX").setdescription("Heat flux magnitude.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("FLUX_X").setdescription("Heat flux in x-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("FLUX_Y").setdescription("Heat flux in y-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("FLUX_Z").setdescription("Heat flux in z-direction.").allowonly([&] () { return _activated && root->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D; }))

	.addoption(ECFOption().setname("GRADIENT").setdescription("Heat gradient magnitude.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("GRADIENT_X").setdescription("Heat gradient in x-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("GRADIENT_Y").setdescription("Heat gradient in y-direction.").allowonly([&] () { return _activated; }))
	.addoption(ECFOption().setname("GRADIENT_Z").setdescription("Heat gradient in z-direction.").allowonly([&] () { return _activated && root->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D; }));
}

HeatTransferOutputSettings::HeatTransferOutputSettings()
{
	basic();

	REGISTER(temperature, ECFMetaData()
			.setdescription({ "Temperature." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
	REGISTER(translation_motions, ECFMetaData()
			.setdescription({ "Translation motions." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
	REGISTER(gradient, ECFMetaData()
			.setdescription({ "Temperature gradient." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));

	REGISTER(flux, ECFMetaData()
			.setdescription({ "Temperature flux." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
	REGISTER(htc, ECFMetaData()
			.setdescription({ "Heat transfer coefficient." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
	REGISTER(phase, ECFMetaData()
			.setdescription({ "Material phase." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
	REGISTER(latent_heat, ECFMetaData()
			.setdescription({ "Latent heat." })
			.setdatatype({ ECFDataType::BOOL })
			.allowonly([&] () { return _activated; }));
}

HeatTransferGlobalSettings::HeatTransferGlobalSettings(ECFObject *ecfdescription)
{
	stabilization = STABILIZATION::SUPG;
	REGISTER(stabilization, ECFMetaData()
			.setdescription({ "Inconsistent stabilization" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SUPG").setdescription("SUPG stabilization"))
			.addoption(ECFOption().setname("CAU").setdescription("CAU stabilization"))
			.setform());

	kernel = KERNEL::OLD;
	REGISTER(kernel, ECFMetaData()
			.setdescription({ "Kernel assembler" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("OLD").setdescription("Old (hopefully) slow assembler."))
			.addoption(ECFOption().setname("OPT").setdescription("Opt assembler.")));

	sigma = 0;
	REGISTER(sigma, ECFMetaData()
			.setdescription({ "Inconsistent stabilization parameter" })
			.setdatatype({ ECFDataType::FLOAT })
			.setform());

	init_temp_respect_bc = true;
	REGISTER(init_temp_respect_bc, ECFMetaData()
			.setdescription({ "Initial temperature follows BC" })
			.setdatatype({ ECFDataType::BOOL })
			.setform());

	diffusion_split = false;
	REGISTER(diffusion_split, ECFMetaData()
			.setdescription({ "Thermal shock stabilization" })
			.setdatatype({ ECFDataType::BOOL })
			.setform());
}

HeatTransferConfiguration::HeatTransferConfiguration(DIMENSION d)
: PhysicsConfiguration(d, MaterialConfiguration::PHYSICAL_MODEL::THERMAL),
  HeatTransferGlobalSettings(ecfdescription),
  dimension(d)
{
	REGISTER(load_steps_settings, ECFMetaData()
			.setdescription({ "Settings for each load step", "LoadStep" })
			.setdatatype({ ECFDataType::LOAD_STEP })
			.setpattern({ "1" })
			.setdynamic(),
			&dimension);
}


