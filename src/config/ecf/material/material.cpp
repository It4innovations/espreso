
#include "material.h"

#include "config/configuration.hpp"

using namespace espreso;

MaterialBaseConfiguration::MaterialBaseConfiguration(DIMENSION *D, PHYSICAL_MODEL physicalModel, bool *phase_change)
: physical_model(physicalModel), material_model(MATERIAL_MODEL::LINEAR_ELASTIC),
  coordinate_system(D),
  density(ECFMetaData::getmaterialvariables()), heat_capacity(ECFMetaData::getmaterialvariables()),
  linear_elastic_properties(D), hyper_elastic_properties(D), thermal_expansion(D),
  thermal_conductivity(D),
  _phase_change(phase_change)
{
	REGISTER(coordinate_system, ECFMetaData()
			.setdescription({ "Coordinate system" })
			.allowonly([&] () { return !*_phase_change; }));

	density.value = heat_capacity.value = "0";
	ecfdescription->registerParameter("dens", density, ECFMetaData()
			.setdescription({ "Density" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit(SIUnit(-3, 1, 0, 0, 0, 0, 0))
			.allowonly([&] () { return !*_phase_change; }));

	ecfdescription->registerParameter("CP", heat_capacity, ECFMetaData()
			.setname("Heat capacity")
			.setdescription({ "Heat capacity" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return !*_phase_change; }));

	if (physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) {
		REGISTER(linear_elastic_properties, ECFMetaData()
				.setdescription({ "Linear elasticity" })
				.allowonly([&] () { return
						(!*_phase_change) &&
						(physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) &&
						(material_model == MATERIAL_MODEL::LINEAR_ELASTIC); }));

		REGISTER(hyper_elastic_properties, ECFMetaData()
				.setdescription({ "Hyper elasticity" })
				.allowonly([&] () { return
						(!*_phase_change) &&
						(physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) &&
						(material_model == MATERIAL_MODEL::HYPER_ELASTIC); }));
		REGISTER(thermal_expansion, ECFMetaData()
				.setdescription({ "Thermal expansion" })
				.allowonly([&] () {
					return (!*_phase_change) && (physical_model & PHYSICAL_MODEL::THERMAL);
				}));
	}

	if (physical_model & PHYSICAL_MODEL::THERMAL) {
		REGISTER(thermal_conductivity, ECFMetaData()
				.setdescription({ "Thermal conductivity" })
				.allowonly([&] () {
					return (!*_phase_change) && (physical_model & PHYSICAL_MODEL::THERMAL);
				}));
	}
}

MaterialConfiguration::MaterialConfiguration(DIMENSION *D, PHYSICAL_MODEL physicalModel)
: MaterialBaseConfiguration(D, physicalModel, &phase_change),
  phase_change(false)
{
	REGISTER(dimension, ECFMetaData()
				.setdescription({"Dimension"})
				.setdatatype({ ECFDataType::OPTION })
				.addoption(ECFOption().setname("D1").setdescription("D1"))
				.addoption(ECFOption().setname("D2").setdescription("D2"))
				.addoption(ECFOption().setname("D3").setdescription("D3"))
				.addoption(ECFOption().setname("Z").setdescription("Z"))
				.allowonly([&] () { return false; })
				.addconstraint(ECFFalseCondition()));

	name = "";
	REGISTER(name, ECFMetaData()
			 .setdescription({ "Name" })
			 .setdatatype( { ECFDataType::STRING } ));

	description = "";
	REGISTER(description, ECFMetaData()
			 .setdescription({ "Description" })
			 .setdatatype( { ECFDataType::STRING } ));

	ecfdescription->addSpace();

	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Physical model" })
			.setdatatype({ ECFDataType::ENUM_FLAGS })
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by HEAT TRANSFER.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::THERMAL; }))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS").setdescription("One of models used by STRUCTURAL MECHANICS.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; })));

	REGISTER(material_model, ECFMetaData()
			.setdescription({ "Material model" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("Linear elastic model.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; }))
			.addoption(ECFOption().setname("HYPER_ELASTIC").setdescription("Hyper elastic model.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; })));

	ecfdescription->addSeparator();

	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);

	REGISTER(phase_change, ECFMetaData()
	.setdescription({ "Phase change" })
	.setdatatype({ ECFDataType::BOOL }));

	smooth_step_order = 1;
	REGISTER(smooth_step_order, ECFMetaData()
			.setdescription({ "Smooth step order" })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.allowonly([&] () { return phase_change; }));

	latent_heat = transition_interval = phase_change_temperature = 0;
	REGISTER(latent_heat, ECFMetaData()
			.setdescription({ "Latent heat" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return phase_change; }));
	REGISTER(transition_interval, ECFMetaData()
				.setdescription({ "Transition interval" })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; }));
	REGISTER(phase_change_temperature, ECFMetaData()
				.setdescription({ "Phase change temperature" })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; }));

	REGISTER(phases, ECFMetaData()
			.setdescription({ "Phase change settings.", "Phase settings" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.setpattern({ "1" })
			.allowonly([&] () { return phase_change; }),
			D, physical_model, &phase_change);

	ecfdescription->getParameter(&phases)->getParameter("1");
	ecfdescription->getParameter(&phases)->getParameter("2");
}
