
#include "material.h"

#include "config/configuration.hpp"
#include "config/conditions.h"

using namespace espreso;

MaterialBaseConfiguration::MaterialBaseConfiguration(DIMENSION *D, PHYSICAL_MODEL physicalModel, bool *phase_change)
: physical_model(physicalModel), material_model(MATERIAL_MODEL::LINEAR_ELASTIC),
  coordinate_system(D), density(ECFExpression::Scope::EGPS), speed_of_sound(ECFExpression::Scope::EGPS), heat_capacity(ECFExpression::Scope::EGPS),
  linear_elastic_properties(D), hyper_elastic_properties(D), thermal_expansion(D),
  thermal_conductivity(D),
  _phase_change(phase_change)
{
	REGISTER(coordinate_system, ECFMetaData()
			.setdescription({ "Coordinate system" })
			.setcollapsed()
			.allowonly([&] () { return !*_phase_change; })
			.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false)));

	density.value = heat_capacity.value = "0";
	ecfdescription->registerParameter("dens", density, ECFMetaData()
			.setdescription({ "Density" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return !*_phase_change; })
			.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false)));

	speed_of_sound.value = "343"; // [m/s] default for air
	ecfdescription->registerParameter("speed_of_sound", speed_of_sound, ECFMetaData()
			.setdescription({ "Speed of sound" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return !*_phase_change && PHYSICAL_MODEL::ACOUSTICS; }));
			/*
			.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false)));
			*/

	ecfdescription->registerParameter("CP", heat_capacity, ECFMetaData()
			.setname("Heat capacity")
			.setdescription({ "Heat capacity" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.setunit(Unit().add(Unit::UnitLibrary::KILOGRAM, 1)
						.add(Unit::UnitLibrary::METRE, 2)
						.add(Unit::UnitLibrary::SECOND, -2)
						.add(Unit::UnitLibrary::KELVIN, -1))
			.allowonly([&] () { return !*_phase_change; })
			.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false)));

	if (physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) {
		REGISTER(linear_elastic_properties, ECFMetaData()
				.setdescription({ "Linear elasticity" })
				.allowonly([&] () { return
						(_phase_change == NULL || !*_phase_change) &&
						(physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) &&
						(material_model == MATERIAL_MODEL::LINEAR_ELASTIC); })
				.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false) & ECFCondition(material_model, ECFCondition::EQUALS, MATERIAL_MODEL::LINEAR_ELASTIC)));

		REGISTER(hyper_elastic_properties, ECFMetaData()
				.setdescription({ "Hyper elasticity" })
				.allowonly([&] () { return
						(_phase_change == NULL || !*_phase_change) &&
						(physical_model & PHYSICAL_MODEL::STRUCTURAL_MECHANICS) &&
						(material_model == MATERIAL_MODEL::HYPER_ELASTIC); })
				.addconstraint(ECFCondition(*_phase_change, ECFCondition::EQUALS, false) & ECFCondition(material_model, ECFCondition::EQUALS, MATERIAL_MODEL::HYPER_ELASTIC)));
		
		REGISTER(thermal_expansion, ECFMetaData()
				.setdescription({ "Thermal expansion" })
				.allowonly([&] () {
					return (_phase_change == NULL || !*_phase_change) && (physical_model & PHYSICAL_MODEL::THERMAL);
				}));
	}

	if (physical_model & PHYSICAL_MODEL::THERMAL) {
		REGISTER(thermal_conductivity, ECFMetaData()
				.setdescription({ "Thermal conductivity" })
				.noexport()
				.allowonly([&] () {
					return (_phase_change == NULL || !*_phase_change) && (physical_model & PHYSICAL_MODEL::THERMAL);
				}));
	}
}

MaterialConfiguration::MaterialConfiguration(DIMENSION *D, PHYSICAL_MODEL physicalModel)
: MaterialBaseConfiguration(D, physicalModel, &phase_change),
  phase_change(false)
{
	name = "";
	REGISTER(name, ECFMetaData()
			 .setdescription({ "Name" })
			 .setdatatype( { ECFDataType::STRING } ));

	description = "";
	REGISTER(description, ECFMetaData()
			 .setdescription({ "Description" })
			 .setdatatype( { ECFDataType::STRING } ));

	ecfdescription->addSpace()->metadata.noexport();

	REGISTER(physical_model, ECFMetaData()
			.setdescription({ "Physical model" })
			.setdatatype({ ECFDataType::ENUM_FLAGS })
			.noexport()
			.addoption(ECFOption().setname("THERMAL").setdescription("Model used by HEAT TRANSFER.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::THERMAL; }))
			.addoption(ECFOption().setname("STRUCTURAL_MECHANICS").setdescription("One of models used by STRUCTURAL MECHANICS.")
					.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; })));
	
	if (physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS)
	{
		REGISTER(material_model, ECFMetaData()
				.setdescription({ "Material model" })
				.setdatatype({ ECFDataType::OPTION })
				.addoption(ECFOption().setname("LINEAR_ELASTIC").setdescription("Linear elastic model.")
						.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; }))
				.addoption(ECFOption().setname("HYPER_ELASTIC").setdescription("Hyper elastic model.")
						.allowonly([&] () { return physicalModel & PHYSICAL_MODEL::STRUCTURAL_MECHANICS; })));
	}
	ecfdescription->addSeparator();

	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);
	// ecfdescription->moveLastBefore(ecfdescription->parameters.front()->name);

	REGISTER(phase_change, ECFMetaData()
	.setdescription({ "Phase change" })
	.setdatatype({ ECFDataType::BOOL }));

	smooth_step_order = 1;
	REGISTER(smooth_step_order, ECFMetaData()
			.setdescription({ "Smooth step order" })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.setrange(0, 100)
			.allowonly([&] () { return phase_change; })
			.addconstraint(ECFCondition(phase_change, ECFCondition::EQUALS, true)));

	latent_heat = transition_interval = phase_change_temperature = 0;
	REGISTER(latent_heat, ECFMetaData()
			.setdescription({ "Latent heat" })
			.setdatatype({ ECFDataType::FLOAT })
			.allowonly([&] () { return phase_change; })
			.addconstraint(ECFCondition(phase_change, ECFCondition::EQUALS, true)));
	REGISTER(transition_interval, ECFMetaData()
				.setdescription({ "Transition interval" })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; })
				.addconstraint(ECFCondition(phase_change, ECFCondition::EQUALS, true)));
	REGISTER(phase_change_temperature, ECFMetaData()
				.setdescription({ "Phase change temperature" })
				.setdatatype({ ECFDataType::FLOAT })
				.allowonly([&] () { return phase_change; })
				.addconstraint(ECFCondition(phase_change, ECFCondition::EQUALS, true)));

	REGISTER(phases, ECFMetaData()
			.setdescription({ "Phase change settings.", "Phase settings" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER })
			.setpattern({ "1" })
			.setgroup()
			.allowonly([&] () { return phase_change; })
			.addconstraint(ECFCondition(phase_change, ECFCondition::EQUALS, true)),
			D, physical_model, &phase_change);

	auto removePhaseConstraints = [&] (ECFParameter* object) {
		ECFObject* obj = static_cast<ECFObject*>(object);
		obj->forEachParameters([&] (ECFParameter* p){
			if (p->metadata.condition->isset())
			{
				if (p->metadata.condition->match(&phase_change))
				{
					p->metadata.removeconstraint();
				}
			}
		}, true, false, true);
	};
	ECFParameter* ph = ecfdescription->getParameter(&phases)->getParameter("1");
	ph->metadata.name = "Phase 1";
	ph->metadata.removeconstraint();
	ph->metadata.setcollapsed();
	removePhaseConstraints(ph);
	ph = ecfdescription->getParameter(&phases)->getParameter("2");
	ph->metadata.name = "Phase 2";
	ph->metadata.removeconstraint();
	ph->metadata.setcollapsed();
	removePhaseConstraints(ph);
}
