
#include "elasticity.mover.h"
#include "heattransfer.mover.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

ElasticityElementIterator::ElasticityElementIterator(ElasticityElementIterator *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration, int dimension, bool omitTemp)
: coordinates(dimension), displacement(dimension, 0), cos(dimension, 0), sin(dimension, 0), acceleration(dimension, 0), thickness(1, 1),
  temperature(1, 273.15), initialTemperature(1, 273.15), angularVelocity(3, 0),
  phase(dimension, 0), displacementAmplitude(dimension, 0), velocity(dimension, 0), velocityAmplitude(dimension, 0), accelerationAmplitude(dimension, 0),
  principalStress(3, 0), componentStress(6, 0), vonMisesStress(1, 0),
  designVariable(1, configuration.topology_optimization_settings.constraint.value), complianceDerivation(1, 0), compliance(1, 0),
  minDesignVariable(configuration.topology_optimization_settings.solver_settings.min_density),
  penaltyFactor(configuration.topology_optimization_settings.solver_settings.penalty_factor),
  harmonic(configuration.type == LoadStepSolverConfiguration::TYPE::HARMONIC), largeDisplacement(configuration.large_displacement)
{
	coordinates.set(info::mesh->nodes->coordinates, info::mesh->elements->procNodes);

	displacement   .setInput(configuration.displacement    , nodeparams  , info::mesh->elements->procNodes);
	temperature    .setInput(configuration.temperature     , kernelparams, info::mesh->elements->procNodes);
	acceleration   .setInput(configuration.acceleration    , kernelparams, info::mesh->elements->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
	angularVelocity.setInput(configuration.angular_velocity, kernelparams, info::mesh->elements->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
	if (harmonic) {
		cos.setInput(configuration.displacement, nodeparams, info::mesh->elements->procNodes);
		sin.setInput(configuration.displacement, nodeparams, info::mesh->elements->procNodes);
	}

	if (previous) {
		displacement.output = previous->displacement.output;
		thickness           = previous->thickness; previous->thickness.kernel.values = NULL;
		initialTemperature  = previous->initialTemperature; previous->initialTemperature.kernel.values = NULL;
		if (harmonic) {
			cos.output = previous->cos.output; previous->cos.kernel.values = NULL;
			sin.output = previous->sin.output; previous->sin.kernel.values = NULL;

			phase                .output = previous->phase.output;
			displacementAmplitude.output = previous->displacementAmplitude.output;
			velocity             .output = previous->velocity.output;
			velocityAmplitude    .output = previous->velocityAmplitude.output;
			acceleration         .output = previous->acceleration.output;
			accelerationAmplitude.output = previous->accelerationAmplitude.output;
		}

		principalStress.output = previous->principalStress.output;
		componentStress.output = previous->componentStress.output;
		vonMisesStress .output = previous->vonMisesStress.output;
	} else {
		displacement.setOutput(NamedData::DataType::VECTOR, "DISPLACEMENT", true, step::TYPE::TIME | step::TYPE::FTT);
		if (harmonic) {
			cos     .setOutput(NamedData::DataType::VECTOR, "DISPLACEMENT_COS", true, step::TYPE::FREQUENCY);
			sin     .setOutput(NamedData::DataType::VECTOR, "DISPLACEMENT_SIN", true, step::TYPE::FREQUENCY);

			phase                .setOutput(NamedData::DataType::SCALAR, "PHASE"                 , info::ecf->output.results_selection.phase, step::TYPE::FREQUENCY);
			displacementAmplitude.setOutput(NamedData::DataType::SCALAR, "DISPLACEMENT_AMPLITUDE", info::ecf->output.results_selection.displacement, step::TYPE::FREQUENCY);
			velocity             .setOutput(NamedData::DataType::VECTOR, "VELOCITY"              , true, step::TYPE::TIME | step::TYPE::FTT);
			velocityAmplitude    .setOutput(NamedData::DataType::SCALAR, "VELOCITY_AMPLITUDE"    , true, step::TYPE::FREQUENCY);
			acceleration         .setOutput(NamedData::DataType::VECTOR, "ACCELERATION"          , true, step::TYPE::TIME | step::TYPE::FTT);
			accelerationAmplitude.setOutput(NamedData::DataType::SCALAR, "ACCELERATION_AMPLITUDE", true, step::TYPE::FREQUENCY);
		}

		if (configuration.topology_optimization) {
			designVariable.setOutput(NamedData::DataType::SCALAR, "DESIGN_VARIABLE");
			complianceDerivation.setOutput(NamedData::DataType::SCALAR, "COMPLIANCE_DERIVATION");
			compliance.setOutput(NamedData::DataType::SCALAR, "COMPLIANCE");
		}

		if (dimension == 2) {
			thickness.setInput(physics.thickness, kernelparams, info::mesh->elements->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
			thickness.setOutput(NamedData::DataType::SCALAR, "THICKNESS", info::ecf->output.results_selection.thickness);
			initialTemperature.setInput(physics.initial_temperature, kernelparams, info::mesh->elements->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
		}

		if (dimension == 3) {
			initialTemperature.setInput(physics.initial_temperature, kernelparams, info::mesh->elements->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
		}

		principalStress.setOutput(NamedData::DataType::NUMBERED   , "PRINCIPAL_STRESS", info::ecf->output.results_selection.stress, step::TYPE::TIME | step::TYPE::FTT);
		componentStress.setOutput(NamedData::DataType::TENSOR_SYMM, "COMPONENT_STRESS", info::ecf->output.results_selection.stress, step::TYPE::TIME | step::TYPE::FTT);
		vonMisesStress .setOutput(NamedData::DataType::SCALAR     , "VON_MISES_STRESS", info::ecf->output.results_selection.stress, step::TYPE::TIME | step::TYPE::FTT);
	}

	nodeparams.coords(3, &info::mesh->nodes->coordinates->datatarray().data()->x);
	kernelparams.coords(dimension, coordinates.kernel.values->datatarray().data());
	kernelparams.temp(temperature.kernel.values->datatarray().data());

	int d = dimension;
	registerSolution(displacement   , { d, 0 });
	if (harmonic) {
		registerSolution(cos,         { d, 0 });
		registerSolution(sin,         { d, 0 });
	}

	registerInput(coordinates       , { d, 0 });
	registerInput(thickness         , { 1, 0 });
	registerInput(initialTemperature, { 1, 0 });
	if (!omitTemp) {
		registerInput(temperature       , { 1, 0 });
	}
	registerInput(acceleration      , { d, 0 });
	registerInput(angularVelocity   , { 3, 0 });

	if (harmonic) {
		registerOutput(phase                , { d, 0 });
		registerOutput(displacementAmplitude, { d, 0 });
		registerOutput(velocity             , { d, 0 });
		registerOutput(velocityAmplitude    , { d, 0 });
		registerOutput(accelerationAmplitude, { d, 0 });
	}

	registerOutput(principalStress     , { 0, 3 });
	registerOutput(componentStress     , { 0, 6 });
	registerOutput(vonMisesStress      , { 0, 1 });
	registerOutput(designVariable      , { 0, 1 });
	registerOutput(complianceDerivation, { 0, 1 });
	registerOutput(compliance          , { 0, 1 });
}

ElasticityElementIterator::ElasticityElementIterator(HeatTransferElementIterator *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration, int dimension)
: ElasticityElementIterator((ElasticityElementIterator*)NULL, physics, gsettings, configuration, dimension, true)
{
	temperature.output = previous->temperature.output;
	registerSolution(temperature, { 1, 0 });
}

ElasticityBoundaryIterator::RotatingForce::RotatingForce(int dimension)
: axis(dimension, 0), radius(1, 0), mass(1, 0), phaseAngle(1, 0), location(1, 0)
{

}

ElasticityBoundaryIterator::ElasticityBoundaryIterator(BoundaryRegionStore *region, ElasticityElementIterator &iterator, StructuralMechanicsLoadStepConfiguration &configuration, int dimension)
: coordinates(dimension), thickness(1, 1), normalPressure(1, 0), force(dimension, 0), harmonicForce(dimension, 0), rotatingForce(dimension)
{
	if (region->dimension == 0) {
		coordinates.set(info::mesh->nodes->coordinates, region->nodes);
		force.setInput(configuration.force, region->name, kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);

		if (configuration.harmonic_force.find(region->name) != configuration.harmonic_force.end()) {
			harmonicForce.setInput(configuration.harmonic_force, region->name, kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
//			harmonicForce.phase    .setInput(configuration.harmonic_force.find(region->name)->second.phase    , kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
		}
		if (configuration.rotating_force.find(region->name) != configuration.rotating_force.end()) {
			rotatingForce.axis      .setInput(configuration.rotating_force.find(region->name)->second.rotation_axis        , kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
			rotatingForce.radius    .setInput(configuration.rotating_force.find(region->name)->second.rotation_radius      , kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
			rotatingForce.mass      .setInput(configuration.rotating_force.find(region->name)->second.unbalance_mass       , kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
			rotatingForce.phaseAngle.setInput(configuration.rotating_force.find(region->name)->second.unbalance_phase_angle, kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
			rotatingForce.location  .setInput(configuration.rotating_force.find(region->name)->second.location             , kernelparams, region->nodes, MoverParameter::Properties::ALLOW_CONSTANT);
		}
	} else {
		coordinates.set(info::mesh->nodes->coordinates, region->procNodes);
		if (dimension == 2) {
			thickness.setInput(iterator.thickness.output, region->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
		}
		normalPressure.setInput(configuration.normal_pressure, region->name, kernelparams, region->procNodes, MoverParameter::Properties::ALLOW_CONSTANT);
	}

	kernelparams.coords(dimension, coordinates.kernel.values->datatarray().data());

	int d = dimension;
	registerInput(coordinates   , { d, 0 });
	registerInput(thickness     , { 1, 0 }); // FIXME: set correct thickness
	registerInput(normalPressure, { 1, 0 });
	registerInput(force         , { d, 0 });

	registerInput(harmonicForce.magnitude , { d, 0 });
	registerInput(harmonicForce.phase     , { d, 0 });
	registerInput(rotatingForce.axis      , { d, 0 });
	registerInput(rotatingForce.radius    , { 1, 0 });
	registerInput(rotatingForce.mass      , { 1, 0 });
	registerInput(rotatingForce.phaseAngle, { 1, 0 });
	registerInput(rotatingForce.location  , { 1, 0 });
}

bool ElasticityBoundaryIterator::hasSettings()
{
	return normalPressure.input.ecf || force.input.ecf || harmonicForce.magnitude.input.ecf || rotatingForce.mass.input.ecf;
}
