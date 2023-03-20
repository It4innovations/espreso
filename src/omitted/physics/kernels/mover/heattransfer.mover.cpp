
#include "heattransfer.mover.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "config/ecf/physics/heattransfer.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

HeatTransferElementIterator::HeatTransferElementIterator(HeatTransferElementIterator *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration, int dimension)
: coordinates(dimension), temperature(1, 273.15), thickness(1, 1), motion(dimension, 0),
  initialTemperature(1, 0), heatSource(1, 0), bioHeat(1, 0), bioHeatDerivation(1, 0),
  phase(1, 0), latentHeat(1, 0), gradient(dimension, 0), flux(dimension, 0), htc(1, 0)
{
	switch (configuration.human_thermoregulation_system.activity_level_unit){
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::STANDING:activity_level_unit=1.1; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::EASY_WORK:activity_level_unit=2.6; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::HARD_WORK:activity_level_unit=4.3; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::WALKING_0_PRCT:activity_level_unit=2.4; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::WALKING_5_PRCT:activity_level_unit=3.5; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::WALKING_15_PRCT:activity_level_unit=5.7; break;
		case HumanThermoregulationSystem::ACTIVITY_LEVEL_UNIT::TEACHER:activity_level_unit=1.5; break;
	}

	coordinates.set(info::mesh->nodes->coordinates, info::mesh->elements->nodes);

	temperature      .setInput(configuration.temperature            , nodeparams  , info::mesh->elements->nodes);
	heatSource       .setInput(configuration.heat_source            , kernelparams, info::mesh->elements->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
	motion           .setInput(configuration.translation_motions    , kernelparams, info::mesh->elements->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
//	bioHeat          .setInput(configuration.bio_heat               , kernelparams, info::mesh->elements->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
//	bioHeatDerivation.setInput(configuration.bio_heat               , kernelparams, info::mesh->elements->nodes, MoverParameter::Properties::ALLOW_CONSTANT);

	if (previous) {
		temperature.output = previous->temperature.output;
		motion.output      = previous->motion.output;

		phase      .output = previous->phase.output;
		latentHeat .output = previous->latentHeat.output;
		gradient   .output = previous->gradient.output;
		flux       .output = previous->flux.output;
		htc        .output = previous->htc.output;

		initialTemperature = previous->initialTemperature; previous->initialTemperature.kernel.values = NULL;
		thickness          = previous->thickness; previous->thickness.kernel.values = NULL;
	} else {
		temperature .setOutput(NamedData::DataType::SCALAR, "TEMPERATURE");
		motion      .setOutput(NamedData::DataType::VECTOR, "TRANSLATION_MOTION", info::ecf->output.results_selection.translation_motions);

		phase       .setOutput(NamedData::DataType::SCALAR, "PHASE"             , info::ecf->output.results_selection.phase && info::mesh->hasPhaseChange());
		latentHeat  .setOutput(NamedData::DataType::SCALAR, "LATENT_HEAT"       , info::ecf->output.results_selection.latent_heat && info::mesh->hasPhaseChange());
		gradient    .setOutput(NamedData::DataType::VECTOR, "GRADIENT"          , info::ecf->output.results_selection.gradient);
		flux        .setOutput(NamedData::DataType::VECTOR, "FLUX"              , info::ecf->output.results_selection.flux);
		htc         .setOutput(NamedData::DataType::SCALAR, "HTC"               , info::ecf->output.results_selection.htc);
		if (dimension == 2) {
			initialTemperature.setInput(physics.initial_temperature, kernelparams, info::mesh->elements->nodes);
			thickness.setInput(physics.thickness, kernelparams, info::mesh->elements->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
			thickness.setOutput(NamedData::DataType::SCALAR, "THICKNESS"        , info::ecf->output.results_selection.thickness);
		}
		if (dimension == 3) {
			initialTemperature.setInput(physics.initial_temperature, kernelparams, info::mesh->elements->nodes);
		}
	}

	for (auto it = configuration.bio_heat.begin(); it != configuration.bio_heat.end(); ++it) {
		info::mesh->eregion(it->first);
	}

	nodeparams.coords(3, &info::mesh->nodes->coordinates->datatarray().data()->x);
	nodeparams.temp(temperature.output.data->data.data());
	kernelparams.coords(dimension, coordinates.kernel.values->datatarray().data());
	kernelparams.temp(temperature.kernel.values->datatarray().data());
	kernelparams.inittemp(initialTemperature.kernel.values->datatarray().data());

	int d = dimension;
	registerSolution(temperature    , { 1, 0 });

	registerInput(coordinates       , { d, 0 });
	registerInput(motion            , { d, 0 });
	registerInput(heatSource        , { 1, 0 });
//	registerInput(bioHeat           , { 1, 0 });
//	registerInput(bioHeatDerivation , { 1, 0 });
	registerInput(thickness         , { 1, 0 });

	registerOutput(phase             , { 0, 1 });
	registerOutput(latentHeat        , { 0, 1 });
	registerOutput(gradient          , { 0, d });
	registerOutput(flux              , { 0, d });
	registerOutput(htc               , { 0, 1 });

	initialTemperature.increment    = { 1, 0 };

	// input temperature is defined on node region, hence,
	// we need to set temperature to output.temperature first in the case of init_respect_bc

	if (previous == NULL) {
		now(initialTemperature.input, initialTemperature.kernel);
		temperature.initFrom(physics.initial_temperature);
		if (gsettings.init_temp_respect_bc) {
			now(temperature.input, temperature.output);
			now(temperature.output, temperature.kernel);
		}
	} else {
		if (configuration.update_initial_temperature) {
			now(temperature.output, initialTemperature.kernel);
		}
	}
}

HeatTransferBoundaryIterator::HeatTransferBoundaryIterator(BoundaryRegionStore *region, HeatTransferElementIterator &iterator, HeatTransferLoadStepConfiguration &configuration, int dimension)
: regionArea(0), radiation(false), convection(false),
  coordinates(dimension), temperature(1, 273.15), thickness(1, 1),
  emissivity(1, 0), extemperature(1, 273.15), heatflow(1, 0), heatflux(1, 0), htc(1, 0)
{
	if (region->dimension == 0) {
		return;
	}

	regionArea = region->area;
	convection = configuration.convection.find(region->name) != configuration.convection.end();
	radiation = configuration.diffuse_radiation.find(region->name) != configuration.diffuse_radiation.end();

	coordinates.set(info::mesh->nodes->coordinates, region->elements);
	temperature.setInput(iterator.temperature.output, region->elements);
	thickness  .setInput(iterator.thickness.output  , region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
	heatflow   .setInput(configuration.heat_flow    , region->name, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
	heatflux   .setInput(configuration.heat_flux    , region->name, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);

	htc.setInput(configuration.convection, region->name, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
	htc.output = iterator.htc.output;
	if (convection) {
		extemperature.setInput(configuration.convection.find(region->name)->second.external_temperature, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
	}
	if (radiation) {
		emissivity.setInput(configuration.diffuse_radiation.find(region->name)->second.emissivity, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
		if (!convection) {
			extemperature.setInput(configuration.diffuse_radiation.find(region->name)->second.external_temperature, kernelparams, region->elements, MoverParameter::Properties::ALLOW_CONSTANT);
		}
	}

	kernelparams.coords(dimension, coordinates.kernel.values->datatarray().data());
	kernelparams.temp(temperature.kernel.values->datatarray().data());

	int d = dimension;
	registerInput(coordinates         , { d, 0 });
	registerInput(temperature         , { 1, 0 });
	registerInput(thickness           , { 1, 0 }); // FIXME: set correct thickness
	registerInput(htc                 , { 1, 0 });
	registerInput(emissivity          , { 1, 0 });
	registerInput(extemperature       , { 1, 0 });
	registerInput(heatflow            , { 1, 0 });
	registerInput(heatflux            , { 1, 0 });
}

bool HeatTransferBoundaryIterator::hasSettings()
{
	return convection || radiation || heatflow.input.ecf || heatflux.input.ecf;
}

static double convectionHTC(const ConvectionConfiguration &convection, esint csize, const double *coordinates, double time, double temp);

Move<InputBioHeat, ElementNodeValues>::Move(const InputBioHeat &from, const ElementNodeValues &htc)
: from(from), bioHeat(htc), moved(false)
{

}

void Move<InputBioHeat, ElementNodeValues>::operator()()
{
	if (moved) { return; }

	std::fill(bioHeat.values->datatarray().begin(), bioHeat.values->datatarray().end(), 0);

//	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			Evaluator::Params params;
			for (esint *e = region->elements->datatarray().begin(t), i = 0; e < region->elements->datatarray().end(t); ++e) {
				auto element = info::mesh->elements->nodes->begin() + *e;
				for (auto n = element->begin(); n != element->end(); ++n, ++i) {
					params.coords(from.params->ncoords(), from.params->coords() + *n * from.params->ncoords()).temp(from.params->temp() + *n);
					double temp = from.params->temp()[*n];
					double arte = it->second.arteriar_blood_temperature.evaluator->eval(params);
					double heat = it->second.blood_specific_heat.evaluator->eval(params);
					double dens = it->second.blood_density.evaluator->eval(params);
					double meta = it->second.metabolic_heat_source.evaluator->eval(params);
					double perf = it->second.blood_perfusion.evaluator->eval(params);
					double ref = it->second.reference_temperature.evaluator->eval(params);
					double mu = it->second.mu.evaluator->eval(params);
					bioHeat.values->datatarray()[i] = perf * (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1)) * (arte - temp) + meta + perf * (pow(2, (temp -ref) / 10) - 1);
//					if (i == 150) {
//						printf("temp=%e, bh=%e\n", temp, bioHeat.values->datatarray()[i]);
//					}
				}
			}
		}
	}
}

Move<InputBioHeatDerivation, ElementNodeValues>::Move(const InputBioHeatDerivation &from, const ElementNodeValues &htc)
: from(from), bioHeatDerivation(htc), moved(false)
{

}

void Move<InputBioHeatDerivation, ElementNodeValues>::operator()()
{
	if (moved) { return; }

	std::fill(bioHeatDerivation.values->datatarray().begin(), bioHeatDerivation.values->datatarray().end(), 0);

//	BLOOD_PERFUSION *( 2.^((TEMPERATURE-REFERENCE_TEMPERATURE)/10)*MU*log(2).*(ARTERIAR_BLOOD_TEMPERATURE-TEMPERATURE)/10 -(BLOOD_DENSITY*BLOOD_SPECIFIC_HEAT +MU*(2.^((TEMPERATURE-REFERENCE_TEMPERATURE)/10)-1))) + BLOOD_PERFUSION *2.^((TEMPERATURE-REFERENCE_TEMPERATURE)/10)*log(2)/10;

//	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (auto it = from.ecf->begin(); it != from.ecf->end(); ++it) {
			ElementsRegionStore *region = info::mesh->eregion(it->first);
			Evaluator::Params params;
			for (esint *e = region->elements->datatarray().begin(t), i = 0; e < region->elements->datatarray().end(t); ++e) {
				auto element = info::mesh->elements->nodes->begin() + *e;
				for (auto n = element->begin(); n != element->end(); ++n, ++i) {
					params.coords(from.params->ncoords(), from.params->coords() + *n * from.params->ncoords()).temp(from.params->temp() + *n);
					double temp = from.params->temp()[*n];
					double arte = it->second.arteriar_blood_temperature.evaluator->eval(params);
					double heat = it->second.blood_specific_heat.evaluator->eval(params);
					double dens = it->second.blood_density.evaluator->eval(params);
					double perf = it->second.blood_perfusion.evaluator->eval(params);
					double ref = it->second.reference_temperature.evaluator->eval(params);
					double mu = it->second.mu.evaluator->eval(params);
//					bioHeat.values->datatarray()[i] = perf * (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1)) * (arte - temp) + meta + perf * (pow(2, (temp -ref) / 10) - 1);
					bioHeatDerivation.values->datatarray()[i] = perf * (pow(2, (temp - ref) / 10) * mu * log(2) * (arte - temp) / 10 - (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1))) + perf * pow(2, (temp - ref) / 10) * log(2) / 10;
//					bioHeatDerivation.values->datatarray()[i] = (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1)) * (arte - temp);
//					bioHeatDerivation.values->datatarray()[i] += 0;
//					bioHeatDerivation.values->datatarray()[i] += 0;
//					bioHeatDerivation.values->datatarray()[i] += 0;
//					bioHeatDerivation.values->datatarray()[i] += 0;
//					bioHeatDerivation.values->datatarray()[i] += pow(2, (temp - ref) / 10) + perf * log(2) * pow(2, (temp - ref) / 10) - 1;
//					if (i == 150) {
//						printf("temp=%f, arte=%f, heat=%f, dens=%f, perf=%f, ref=%f, mu=%f\n", temp, arte, heat, dens, perf, ref, mu);
//						printf("AA = pow(2, (temp - ref) / 10) * mu * log10(2) * (arte - temp) / 10 == %f\n", pow(2, (temp - ref) / 10) * mu * log10(2) * (arte - temp) / 10);
//						printf("BB = (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1)) == %f\n", (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1)));
//						printf("perf * (AA - BB) == %f\n", perf * (pow(2, (temp - ref) / 10) * mu * log10(2) * (arte - temp) / 10 - (dens * heat + mu * (pow(2, (temp - ref) / 10) - 1))));
//						printf("perf * pow(2, (temp - ref) / 10) * log10(2) / 10 == %f\n", perf * pow(2, (temp - ref) / 10) * log10(2) / 10);
//						printf("temp=%e, der: %e\n", temp, bioHeatDerivation.values->datatarray()[i]);
//					}
				}
			}
		}
	}
}

Move<InputConvection, ElementNodeValues>::Move(const InputConvection &from, const ElementNodeValues &htc)
: from(from), htc(htc), moved(false)
{

}

void Move<InputConvection, ElementNodeValues>::operator()()
{
	if (from.ecf == NULL) {
		return;
	}

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (size_t i = htc.values->datatarray().distribution()[t]; i < htc.values->datatarray().distribution()[t + 1]; ++i) {
			htc.values->datatarray()[i] = convectionHTC(
					*from.ecf,
					from.params->ncoords(), from.params->coords() + i * from.params->ncoords(),
					step::time.current,
					from.params->temp()[i]);
		}
	}
}

static void convectionMaterialParameters(
		const ConvectionConfiguration &convection,
		esint csize, const double *coordinates, double time, double temp, double extTemp,
		double &rho, double &dynamicViscosity, double &dynamicViscosityTemp, double &heatCapacity, double &thermalConductivity)
{
	double  gas_constant, R, M;

	switch (convection.fluid) {
	case ConvectionConfiguration::FLUID::AIR: {
		gas_constant = 286.9;
		rho = convection.absolute_pressure.evaluator->eval(Evaluator::Params().coords(csize, coordinates).temp(&extTemp).time(time));
		rho /= (gas_constant * extTemp);

		if ((extTemp >=200) && (extTemp <= 1600)){
			heatCapacity = 1047.63657-0.372589265*extTemp+9.45304214E-4*pow(extTemp,2.0)-6.02409443E-7*pow(extTemp,3.0)+1.2858961E-10*pow(extTemp,4.0);
		}else if (extTemp < 200){
			heatCapacity = 1047.63657-0.372589265*200.0+9.45304214E-4*pow(200.0,2.0)-6.02409443E-7*pow(200.0,3.0)+1.2858961E-10*pow(200.0,4.0);
		}else if (extTemp > 1600){
			heatCapacity = 1047.63657-0.372589265*1600.0+9.45304214E-4*pow(1600.0,2.0)-6.02409443E-7*pow(1600.0,3.0)+1.2858961E-10*pow(1600.0,4.0);
		}

		if ((extTemp >=200) && (extTemp <= 1600)){
			thermalConductivity = -0.00227583562+1.15480022E-4*extTemp-7.90252856E-8*pow(extTemp,2.0)+4.11702505E-11*pow(extTemp,3.0)-7.43864331E-15*pow(extTemp,4.0);
		}else if (extTemp < 200){
			thermalConductivity = -0.00227583562+1.15480022E-4*200.0-7.90252856E-8*pow(200.0,2.0)+4.11702505E-11*pow(200.0,3.0)-7.43864331E-15*pow(200.0,4.0);
		}else if (extTemp > 1600){
			thermalConductivity =  -0.00227583562+1.15480022E-4*1600.0-7.90252856E-8*pow(1600.0,2.0)+4.11702505E-11*pow(1600.0,3.0)-7.43864331E-15*pow(1600.0,4.0);
		}

		if ((extTemp >=200) && (extTemp <= 1600)){
		dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * extTemp - 7.69429583E-11 * pow(extTemp,2.0) + 4.6437266E-14 * pow(extTemp,3.0) - 1.06585607E-17 * pow(extTemp,4.0);
		}else if (extTemp < 200){
			dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (extTemp > 1600){
			dynamicViscosity = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}


		if ((temp >=200) && (temp <= 1600)){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * temp - 7.69429583E-11 * pow(temp,2.0) + 4.6437266E-14 * pow(temp,3.0) - 1.06585607E-17 * pow(temp,4.0);
		}else if (temp < 200){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (temp > 1600){
			dynamicViscosityTemp = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}

	}break;

	case ConvectionConfiguration::FLUID::STEAM:{
		R = 8314.4598;
		M = 18.02;
		gas_constant = R/M;

		double pressure = convection.absolute_pressure.evaluator->eval(Evaluator::Params().coords(csize, coordinates).temp(&extTemp).time(time));
		rho = pressure / (gas_constant * extTemp);

		if (extTemp >=100){
			heatCapacity = gas_constant * ( 3.38684 + 34.7498e-4 * extTemp - 6.3547E-06 * pow(extTemp,2.0) + 6.96858E-9 * pow(extTemp,3.0) - 2.50659E-12 * pow(extTemp,4.0) );
		}else if (extTemp < 100){
			heatCapacity = 2044;
		}

		if (extTemp >=100 ){
			thermalConductivity = 5.0225E-8 * pow(extTemp,2.0) + 4.724E-05 *extTemp + 1.2445E-4;
		}else if (extTemp < 100){
			thermalConductivity = 0.0248;
		}

		if (extTemp >=100 ){
			dynamicViscosity = (4.0718E-2 * extTemp - 2.9895)* 1E-6;
		}else if (extTemp < 100){
			dynamicViscosity = 12E-6;
		}
	}break;

	case ConvectionConfiguration::FLUID::WATER:{

		if ((extTemp >=273) && (extTemp <= 283)){
			rho = 972.7584 + 0.2084 *extTemp - 4.0E-4 * pow(extTemp,2.0);
		}else if((extTemp >283) && (extTemp <= 373)){
			rho = 345.28 + 5.749816 * extTemp - 0.0157244 * pow(extTemp,2.0) + 1.264375E-5 * pow(extTemp,3.0);
		}else if (extTemp < 273){
			rho = 972.7584 + 0.2084 *273.0 - 4.0E-4 * pow(273.0,2.0);
		}else if (extTemp > 373){
			rho = 345.28 + 5.749816 * 373.0 - 0.0157244 * pow(373.0,2.0) + 1.264375E-5 * pow(373.0,3.0);
		}


		if ((extTemp >= 265) && (extTemp <= 293)){
			dynamicViscosity = 5.948859 - 0.08236196 * extTemp + 4.287142E-4 * pow(extTemp,2.0) - 9.938045E-7 * pow(extTemp,3.0) + 8.65316E-10 * pow(extTemp,4.0);
		}else if((extTemp >293) && (extTemp <= 353)){
			dynamicViscosity = 	0.410191 - 0.004753985 * extTemp + 2.079795E-5 * pow(extTemp,2.0) - 4.061698E-8 *  pow(extTemp,3.0) + 2.983925E-11 * pow(extTemp,4.0);
		}else if((extTemp >353) && (extTemp <= 423)){
			dynamicViscosity = 0.03625638 - 3.265463E-4 * extTemp + 1.127139E-6 * pow(extTemp,2.0) - 1.75363E-9 * pow(extTemp,3.0) + 1.033976E-12 * pow(extTemp,4.0);
		}else if (extTemp < 265){
			dynamicViscosity = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (extTemp > 423){
			dynamicViscosity = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}


		if ((extTemp >= 275) && (extTemp <= 370)){
			thermalConductivity = -0.9003748 + 0.008387698 * extTemp - 1.118205E-5 * pow(extTemp,2.0);
		}else if (extTemp < 275){
			thermalConductivity = -0.9003748 + 0.008387698 * 275.0 - 1.118205E-5 * pow(275.0,2.0);
		}else if (extTemp > 370){
			thermalConductivity = -0.9003748 + 0.008387698 * 370.0 - 1.118205E-5 * pow(370.0,2.0);
		}

		if ((extTemp >= 293) && (extTemp <= 373)){
			heatCapacity = 4035.841 + 0.492312 * extTemp;
		}else if (extTemp < 293){
			heatCapacity = 4035.841 + 0.492312 * 293.0;
		}else if (extTemp > 373){
			heatCapacity = 4035.841 + 0.492312 * 373.0;
		}


		if ((temp >= 265) && (temp <= 293)){
			dynamicViscosityTemp = 5.948859 - 0.08236196 * temp + 4.287142E-4 * pow(temp,2.0) - 9.938045E-7 * pow(temp,3.0) + 8.65316E-10 * pow(temp,4.0);
		}else if((temp >293) && (temp <= 353)){
			dynamicViscosityTemp = 	0.410191 - 0.004753985 * temp + 2.079795E-5 * pow(temp,2.0) - 4.061698E-8 *  pow(temp,3.0) + 2.983925E-11 * pow(temp,4.0);
		}else if((temp >353) && (temp <= 423)){
			dynamicViscosityTemp = 0.03625638 - 3.265463E-4 * temp + 1.127139E-6 * pow(temp,2.0) - 1.75363E-9 * pow(temp,3.0) + 1.033976E-12 * pow(temp,4.0);
		}else if (temp < 265){
			dynamicViscosityTemp = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (temp > 423){
			dynamicViscosityTemp = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}

	}break;
	case ConvectionConfiguration::FLUID::ENGINE_OIL:{


		if ((extTemp >=273) && (extTemp <= 433)){
			rho = 1068.70404 - 0.6393421 * extTemp + 7.34307359E-5 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			rho = 1068.70404 - 0.6393421 * 273.0 + 7.34307359E-5 * pow(273.0,2.0);
		}else if (extTemp > 433){
			rho = 1068.70404 - 0.6393421 * 433.0 + 7.34307359E-5 * pow(433.0,2.0);
		}

		if ((extTemp >= 273) && (extTemp <= 433)){
			thermalConductivity = 0.192223542 - 2.0637987E-4 * extTemp + 1.54220779E-7 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			thermalConductivity = 0.192223542 - 2.0637987E-4 * 273.0 + 1.54220779E-7 * pow(273.0,2.0);
		}else if (extTemp > 433){
			thermalConductivity =0.192223542 - 2.0637987E-4 * 433.0 + 1.54220779E-7 * pow(433.0,2.0);
		}


		if ((extTemp >= 273) && (extTemp <= 433)){
			heatCapacity = 761.405625 + 3.47685606 * extTemp + 0.00115530303 * pow(extTemp,2.0);
		}else if (extTemp < 273){
			heatCapacity = 761.405625 + 3.47685606 * 273.0 + 0.00115530303 * pow(273.0,2.0);
		}else if (extTemp > 433){
			heatCapacity = 761.405625 + 3.47685606 * 433.0 + 0.00115530303 * pow(433.0,2.0);
		}


		if ((extTemp >= 273) && (extTemp <= 353)){
			dynamicViscosity = 42669.28688622 - 741.1718801282 * extTemp + 5.360521287088 * pow(extTemp,2.0) - 0.02066027676164 * pow(extTemp,3.0) + 4.47491538052E-5 * pow(extTemp,4.0) - 5.164053479202E-8 * pow(extTemp,5.0) + 2.48033770504E-11 * pow(extTemp,6.0);
		}else if ((extTemp > 353) && (extTemp <= 433 )){
			dynamicViscosity = 4.94593941 - 0.0351869631 * extTemp + 8.37935977E-5 * pow(extTemp,2.0) - 6.67125E-8 * pow(extTemp,3.0);

		}else if (extTemp < 273){
			dynamicViscosity = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (extTemp > 433){
			dynamicViscosity = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}

		if ((temp >= 273) && (temp <= 353)){
			dynamicViscosityTemp = 42669.28688622 - 741.1718801282 * temp + 5.360521287088 * pow(temp,2.0) - 0.02066027676164 * pow(temp,3.0) + 4.47491538052E-5 * pow(temp,4.0) - 5.164053479202E-8 * pow(temp,5.0) + 2.48033770504E-11 * pow(temp,6.0);
		}else if ((temp > 353) && (temp <= 433 )){
			dynamicViscosityTemp = 4.94593941 - 0.0351869631 * temp + 8.37935977E-5 * pow(temp,2.0) - 6.67125E-8 * pow(temp,3.0);

		}else if (temp < 273){
			dynamicViscosityTemp = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (temp > 433){
			dynamicViscosityTemp = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}


	}break;
	case ConvectionConfiguration::FLUID::TRANSFORMER_OIL:{

		if ((extTemp >=223) && (extTemp <= 373)){
			rho = 1055.04607 - 0.581753034 * extTemp - 6.40531689E-5 * pow(extTemp,2.0);
		}else if (extTemp < 223){
			rho =  1055.04607 - 0.581753034 * 223.0 - 6.40531689E-5 * pow(223.0,2.0);
		}else if (extTemp > 373){
			rho = 1055.04607 - 0.581753034 * 373.0 - 6.40531689E-5 * pow(373.0,2.0);
		}

		if ((extTemp >= 273) && (extTemp <= 433)){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * extTemp;
		}else if (extTemp < 273){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * 223.0;
		}else if (extTemp > 433){
			thermalConductivity = 0.134299084 - 8.04973822E-5 * 373.0;
		}


		if ((extTemp >= 223) && (extTemp <= 293)){
			heatCapacity = -117056.38 + 1816.76208 * extTemp - 10.305786 * pow(extTemp,2.0) + 0.0256691919 * pow(extTemp,3.0) - 2.36742424E-5 * pow(extTemp,4.0);
		}else if ((extTemp > 293) && (extTemp <= 373 )){
			heatCapacity = -13408.1491 + 123.044152 * extTemp - 0.335401786 * pow(extTemp,2.0) + 3.125E-4 * pow(extTemp,3.0);
		}else if (extTemp < 223){
			heatCapacity = -117056.38 + 1816.76208 * 223.0 - 10.305786 * pow(223.0,2.0) + 0.0256691919 * pow(223.0,3.0) - 2.36742424E-5 * pow(223.0,4.0);
		}else if (extTemp > 373){
			heatCapacity = -13408.1491 + 123.044152 * 373.0 - 0.335401786 * pow(373.0,2.0) + 3.125E-4 * pow(373.0,3.0);
		}


		if ((extTemp >= 243) && (extTemp <= 273)){
			dynamicViscosity = 4492.20229 - 64.7408879 * extTemp + 0.349900959 * pow(extTemp,2.0) - 8.40477E-4 * pow(extTemp,3.0) + 7.57041667E-7 * pow(extTemp,4.0);
		}else if ((extTemp > 273) && (extTemp <= 373 )){
			dynamicViscosity = 91.4524999 - 1.33227058 * extTemp + 0.00777680216 * pow(extTemp,2.0) - 2.27271368E-5 *  pow(extTemp,3.0) + 3.32419673E-8 * pow(extTemp,4.0) - 1.94631023E-11 * pow(extTemp,5.0);
		}else if (extTemp < 243){
			dynamicViscosity = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (extTemp > 373){
			dynamicViscosity = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

		if ((temp >= 243) && (temp <= 273)){
			dynamicViscosityTemp = 4492.20229 - 64.7408879 * temp + 0.349900959 * pow(temp,2.0) - 8.40477E-4 * pow(temp,3.0) + 7.57041667E-7 * pow(temp,4.0);
		}else if ((temp > 273) && (temp <= 373 )){
			dynamicViscosityTemp = 91.4524999 - 1.33227058 * temp + 0.00777680216 * pow(temp,2.0) - 2.27271368E-5 *  pow(temp,3.0) + 3.32419673E-8 * pow(temp,4.0) - 1.94631023E-11 * pow(temp,5.0);
		}else if (temp < 243){
			dynamicViscosityTemp = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (temp > 373){
			dynamicViscosityTemp = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

	}break;
	default:
		eslog::error("Invalid convection fluid type.");
	}
}

static double convectionHTC(
		const ConvectionConfiguration &convection,
		esint csize, const double *coordinates, double time, double temp)
{
	double htc = 0;
	Evaluator::Params params;
	params.coords(csize, coordinates);
	params.temp(&temp);
	params.time(time);
	switch (convection.type) {
	case ConvectionConfiguration::TYPE::USER:
		return convection.heat_transfer_coefficient.evaluator->eval(params);

	case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL: {
		double avgTemp, g, rho, dynamicViscosity, heatCapacity, thermalConductivity, dynamicViscosityTemp, extTemp;

		extTemp = convection.external_temperature.evaluator->eval(params);
		avgTemp = (extTemp + temp) / 2;
		g = 9.81;

		convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);

		switch (convection.variant) {
		case ConvectionConfiguration::VARIANT::INCLINED_WALL: {
			double wallHeight, tiltAngle;
			wallHeight = convection.wall_height.evaluator->eval(params);
			double RaL = pow(rho,2) * g * (1 / avgTemp) * heatCapacity * std::fabs(temp - extTemp  ) * pow(wallHeight, 3.0) / ( thermalConductivity * dynamicViscosity);
			tiltAngle = convection.tilt_angle.evaluator->eval(params);
			tiltAngle *= M_PI / 180.0;
			if (RaL <= 1e9) {
				htc = (thermalConductivity / wallHeight) * (0.68 + (0.67 * cos(tiltAngle) * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermalConductivity / wallHeight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::VERTICAL_WALL: {
			double wallHeight = convection.wall_height.evaluator->eval(params);
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(wallHeight, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (RaL <= 1e9) {
				htc = (thermalConductivity / wallHeight) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermalConductivity / wallHeight) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_UP: {
			double length;
			length = convection.length.evaluator->eval(params);
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(length, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (temp > extTemp) {
				if (RaL <= 1e7) {
					htc = thermalConductivity / length * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermalConductivity / length * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermalConductivity / length * 0.27 * pow(RaL,0.25);
			}

		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_DOWN: {
			double length = convection.length.evaluator->eval(params);
			double RaL = pow(rho,2)	* g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) *pow(length, 3.0)/ ( thermalConductivity * dynamicViscosity);

			if (temp <= extTemp) {
				if (RaL <= 1e7) {
					htc = thermalConductivity / length * 0.54 * pow(RaL,0.25);
				} else {
					htc = thermalConductivity / length * 0.15 * pow(RaL,1.0/3.0);
				}
			} else {
				htc = thermalConductivity / length * 0.27 * pow(RaL,0.25);
			}
		} break;
		case ConvectionConfiguration::VARIANT::HORIZONTAL_CYLINDER:{
			double diameter = convection.diameter.evaluator->eval(params);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
			double Pr = dynamicViscosity * heatCapacity / thermalConductivity;

			if ( RaD > 10e12 ){
				// warning!!!!
				eslog::error("Validated only for RaD <= 10e12.");
			}
			htc = thermalConductivity / diameter * pow( 0.6 + ( 0.387*pow(RaD,1.0/6.0)/ pow( 1 + pow( 0.559/Pr, 9.0/16.0), 8.0/27.0) ) ,2.0);

		} break;
		case ConvectionConfiguration::VARIANT::SPHERE: {
			double diameter = convection.diameter.evaluator->eval(params);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs( temp - extTemp) * pow(diameter, 3.0) / ( thermalConductivity * dynamicViscosity);
			double Pr = dynamicViscosity * heatCapacity / thermalConductivity;

			if ( RaD > 10e11 || Pr < 0.7 ){
				// warning!!!!
				eslog::error("Validated only for RaD <= 10e11 and Pr >= 0.7.");
			}
			htc = thermalConductivity / diameter * pow(2.0 + ( 0.589*pow(RaD,0.25)/ pow( 1 + pow( 0.469/Pr, 9.0/16.0), 4.0/9.0) ) ,2.0);
		} break;
		default:
			eslog::error("Invalid convection variant for EXTERNAL_NATURAL.");
		}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:{

		double avgTemp, g, rho, dynamicViscosity, heatCapacity, thermalConductivity,dynamicViscosityTemp, extTemp;

		extTemp = convection.length.evaluator->eval(params);
		avgTemp = (extTemp + temp) / 2.0;
		g = 9.81;

		convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity );

		switch (convection.variant) {
		case ConvectionConfiguration::VARIANT::PARALLEL_PLATES: {
			double wallHeight, length;
			wallHeight = convection.wall_height.evaluator->eval(params);
			length = convection.length.evaluator->eval(params);
			double H_L = wallHeight / length;
			double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(length,3.0)/ ( thermalConductivity * dynamicViscosity);

			if (( RaL < H_L ) && (temp >  extTemp)) {
				htc = thermalConductivity / wallHeight * ( 1.0 / 24.0 ) * RaL;
			} else {
				double RaL = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs( temp - extTemp) * pow(length,3.0)/ ( thermalConductivity * dynamicViscosity);
				if (RaL <= 1e9) {
					htc = (thermalConductivity / length) * (0.68 + (0.67 * pow(RaL,0.25)) / (pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermalConductivity / length) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		case ConvectionConfiguration::VARIANT::CIRCULAR_TUBE: {
			double diameter, wallHeight;
			diameter = convection.diameter.evaluator->eval(params);
			wallHeight = convection.wall_height.evaluator->eval(params);
			double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
			double H_D = wallHeight / diameter;

			if ( RaD < H_D ){
				htc = thermalConductivity / wallHeight * ( 1.0 / 128.0 ) * RaD;
			}else{

				double RaD = pow(rho,2) * g * (1/avgTemp) * heatCapacity * std::fabs(temp - extTemp) * pow(diameter, 3.0)/ ( thermalConductivity * dynamicViscosity);
				if (RaD <= 1e9) {
					htc = (thermalConductivity / diameter) * (0.68 + (0.67 * pow(RaD,0.25))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermalConductivity / diameter) * pow(0.825 + (0.387 * pow(RaD,1.0/6.0))/(pow( 1+ pow((0.492 * thermalConductivity)/(dynamicViscosity * heatCapacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		} break;
		default:
			eslog::error("Invalid convection variant for INTERNAL_NATURAL.\n");
		}
	} break;

	case ConvectionConfiguration::TYPE::EXTERNAL_FORCED: {

			switch (convection.variant) {
			case ConvectionConfiguration::VARIANT::AVERAGE_PLATE: {

				double avgTemp, rho, dynamicViscosity, heatCapacity, thermalConductivity,dynamicViscosityTemp, extTemp, length, fluidVelocity;
				extTemp = convection.external_temperature.evaluator->eval(params);
				length = convection.length.evaluator->eval(params);
				fluidVelocity = convection.fluid_velocity.evaluator->eval(params);

				avgTemp = (extTemp + temp) / 2.0;

				convectionMaterialParameters(convection, csize, coordinates, time, temp, avgTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);


				double Re = rho * fluidVelocity * length / dynamicViscosity;
				double Pr = dynamicViscosity * heatCapacity / thermalConductivity;
				if (Re <= 5e5) {
					htc = 2 * (thermalConductivity / length) * ((0.3387 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.5)) / (pow(1 + pow(0.0468 / Pr, 2.0 / 3.0), 0.25)));
				} else {
					htc = 2 * (thermalConductivity / length) * pow(Pr, 1.0 / 3.0)	* (0.037 * pow(Re, 0.8) - 871);
				}

			} break;
			default:
				eslog::error("Invalid convection variant for EXTERNAL_FORCED.\n");
			}
	} break;

	case ConvectionConfiguration::TYPE::INTERNAL_FORCED:{

			switch (convection.variant) {
			case ConvectionConfiguration::VARIANT::TUBE: {

				double extTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity, fluidVelocity, diameter;

				extTemp = convection.external_temperature.evaluator->eval(params);
				fluidVelocity = convection.fluid_velocity.evaluator->eval(params);
				diameter = convection.diameter.evaluator->eval(params);

				convectionMaterialParameters(convection, csize, coordinates, time, temp, extTemp, rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity);

				double Re = rho * fluidVelocity * diameter / dynamicViscosity;
				double Pr = dynamicViscosity * heatCapacity / thermalConductivity;
				double n = temp < extTemp ? 0.3 : 0.4;
				htc = thermalConductivity / diameter;
				if (Re <= 2500) {
					htc *= 3.66;
				} else {
					htc *= 0.027 * pow(Re, .8) * pow(Pr, n)	* pow(dynamicViscosity / dynamicViscosityTemp, 0.14);
				}
			}break;

			case ConvectionConfiguration::VARIANT::QUENCH: {

				double text, press, VFRAC, C, g, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity;
				double rhoAIR, dynamic_viscosityAIR, dynamic_viscosity_TAIR, heat_capacityAIR, thermal_conductivityAIR;
				double rhoWATER, dynamic_viscosityWATER, dynamic_viscosity_TWATER, heat_capacityWATER, thermal_conductivityWATER;

				press = convection.absolute_pressure.evaluator->eval(params);
				text = convection.external_temperature.evaluator->eval(params);
				C = convection.experimental_constant.evaluator->eval(params);
				T_AVG = (text + temp) / 2;
				g = 9.81;

				convectionMaterialParameters(convection, csize, coordinates, time, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity);
				htc = C * 0.424 * pow((pow(thermal_conductivity,3.0) * rho * g * (958.35 - rho) * (2257600 + 0.4 * heat_capacity * (temp - (27.952 * log(press) - 222.5304))))/( dynamic_viscosity* (temp - ( 27.952 * log(press) - 222.5304 ) ) * pow( 0.06/(g*( 958.35 - rho )) ,0.5)) ,0.25);
				VFRAC = convection.volume_fraction.evaluator->eval(params);

				ConvectionConfiguration convectionAIR;
				convectionAIR.fluid = ConvectionConfiguration::FLUID::AIR;
				convectionAIR.absolute_pressure = convection.absolute_pressure;

				convectionMaterialParameters(convectionAIR, csize, coordinates, time, temp, T_AVG, rhoAIR, dynamic_viscosityAIR, dynamic_viscosity_TAIR, heat_capacityAIR, thermal_conductivityAIR);


				ConvectionConfiguration convectionWATER;
				convectionAIR.fluid = ConvectionConfiguration::FLUID::WATER;
				convectionAIR.absolute_pressure = convection.absolute_pressure;

				convectionMaterialParameters(convectionAIR, csize, coordinates, time, temp, T_AVG, rhoWATER, dynamic_viscosityWATER, dynamic_viscosity_TWATER, heat_capacityWATER, thermal_conductivityWATER);


				dynamic_viscosity = dynamic_viscosityWATER * VFRAC + dynamic_viscosityAIR*(1-VFRAC);
				heat_capacity = heat_capacityWATER * VFRAC + heat_capacityAIR*(1-VFRAC);
				thermal_conductivity = thermal_conductivityWATER * VFRAC + thermal_conductivityAIR*(1-VFRAC);
				rho = rhoWATER * VFRAC + rhoAIR*(1-VFRAC);

				double length = convection.length.evaluator->eval(params);

				double fvel = convection.fluid_velocity.evaluator->eval(params);
				double Re = rho * fvel * length / dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

 		    	htc =  VFRAC * htc + (1-VFRAC) * (thermal_conductivity / length) * pow(Pr, 1.0 / 3.0)	* (0.228 * pow(Re, 0.731));

			}break;

			case ConvectionConfiguration::VARIANT::QUENCH_PARALLEL: {

				double text, press, VFRAC, C, g, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity;
				double rhoAIR, dynamic_viscosityAIR, dynamic_viscosity_TAIR, heat_capacityAIR, thermal_conductivityAIR;
				double rhoWATER, dynamic_viscosityWATER, dynamic_viscosity_TWATER, heat_capacityWATER, thermal_conductivityWATER;

				press = convection.absolute_pressure.evaluator->eval(params);
				text = convection.external_temperature.evaluator->eval(params);
				C = convection.experimental_constant.evaluator->eval(params);
				T_AVG = (text + temp) / 2;
				g = 9.81;

				convectionMaterialParameters(convection, csize, coordinates, time, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity);
				htc = C * 0.424 * pow((pow(thermal_conductivity,3.0) * rho * g * (958.35 - rho) * (2257600 + 0.4 * heat_capacity * (temp - (27.952 * log(press) - 222.5304))))/( dynamic_viscosity* (temp - ( 27.952 * log(press) - 222.5304 ) ) * pow( 0.06/(g*( 958.35 - rho )) ,0.5)) ,0.25);
				VFRAC = convection.volume_fraction.evaluator->eval(params);

				ConvectionConfiguration convectionAIR;
				convectionAIR.fluid = ConvectionConfiguration::FLUID::AIR;
				convectionAIR.absolute_pressure = convection.absolute_pressure;

				convectionMaterialParameters(convectionAIR, csize, coordinates, time, temp, T_AVG, rhoAIR, dynamic_viscosityAIR, dynamic_viscosity_TAIR, heat_capacityAIR, thermal_conductivityAIR);


				ConvectionConfiguration convectionWATER;
				convectionAIR.fluid = ConvectionConfiguration::FLUID::WATER;
				convectionAIR.absolute_pressure = convection.absolute_pressure;

				convectionMaterialParameters(convectionAIR, csize, coordinates, time, temp, T_AVG, rhoWATER, dynamic_viscosityWATER, dynamic_viscosity_TWATER, heat_capacityWATER, thermal_conductivityWATER);


				dynamic_viscosity = dynamic_viscosityWATER * VFRAC + dynamic_viscosityAIR*(1-VFRAC);
				heat_capacity = heat_capacityWATER * VFRAC + heat_capacityAIR*(1-VFRAC);
				thermal_conductivity = thermal_conductivityWATER * VFRAC + thermal_conductivityAIR*(1-VFRAC);
				rho = rhoWATER * VFRAC + rhoAIR*(1-VFRAC);

				double length = convection.length.evaluator->eval(params);

				double fvel = convection.fluid_velocity.evaluator->eval(params);
				double Re = rho * fvel * length / dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;
				double htcA;

				if (Re <= 5e5) {
					htcA = 2.0 * (thermal_conductivity / length) * ((0.3387 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.5)) / (pow( 1 + 0.0468 / pow( Pr, 2.0 / 3.0), 0.25)));
				} else {
					htcA = 2.0 * (thermal_conductivity / length) * pow(Pr, 1.0 / 3.0)	* (0.037 * pow(Re, 0.8) - 871);
				}

 		    	htc =  VFRAC * htc + (1-VFRAC) * htcA + 5.0;

			}break;



			default:
				eslog::error("Invalid convection variant for EXTERNAL_FORCED.\n");
			}
	}break;

	default:
		eslog::error("Invalid convection TYPE.\n");
	}

	return htc;
}

