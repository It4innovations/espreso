
#include "structuralmechanics.h"
#include "assembler.hpp"

#include "basis/expression/variable.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/assembler/operators/operators.h"
#include "analysis/scheme/steadystate.h"

#include <numeric>
#include <algorithm>

#include "basis/utilities/print.h"

using namespace espreso;

NodeData* StructuralMechanics::Results::displacement = nullptr;

StructuralMechanics::StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{

}

void StructuralMechanics::initParameters()
{
	if (Results::displacement == nullptr) {
		Results::displacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT");
		Variable::list.node["DISPLACEMENT_X"] = new OutputVariable(Results::displacement, 0, info::mesh->dimension);
		Variable::list.node["DISPLACEMENT_Y"] = new OutputVariable(Results::displacement, 1, info::mesh->dimension);
		if (info::mesh->dimension == 3) {
			Variable::list.node["DISPLACEMENT_Z"] = new OutputVariable(Results::displacement, 2, info::mesh->dimension);
		}
	}
}

bool StructuralMechanics::initTemperature()
{
	// This code has to solve problem that initial temperature is set to elements regions, but we need it in nodes
	// 1. settings -> nodeInitialTemperature
	// 2. nodeInitialTemperature -> initialTemperature (here we average values)
	// 3. Dirichlet -> initialTemperature (if 'init_temp_respect_bc')
	// 4. initialTemperature -> nodeInitialTemperature (TODO: check the correction with TB)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool correct = true;

//	if (configuration.temperature.size()) {
//		correct &= examineBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature, temperature.node.externalValues);
//		fromExpression(*this, temperature.node, temperature.node.externalValues);
//	}
//
//	if (settings.init_temp_respect_bc) {
//		temp.initial.node.setConstness(false);
//	}
//	examineElementParameter("INITIAL TEMPERATURE", settings.initial_temperature, temp.initial.node.externalValues);
//	fromExpression(*this, temp.initial.node, temp.initial.node.externalValues);
//	_evaluate();
//
////	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
////		StructuralMechanics::insertParameters(it->second.evaluator);
////	}
////
////	temp.initial.node.builder->buildAndExecute(*this);
//
//	averageEnodesToNodes(temp.initial.node, *ParametersTemperature::Initial::output);
//	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
//
////	if (info::mesh->dimension == 2 && info::ecf->heat_transfer_2d.init_temp_respect_bc) {
////		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::Initial::output, "SET INITIAL TEMPERATURE ACCORDING TO DIRICHLET").buildAndExecute(*this);
////	}
////	if (info::mesh->dimension == 3 && info::ecf->heat_transfer_3d.init_temp_respect_bc) {
////		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::Initial::output, "SET INITIAL TEMPERATURE ACCORDING TO DIRICHLET").buildAndExecute(*this);
////	}
////	CopyNodesToElementsNodes(*ParametersTemperature::Initial::output, temp.initial.node, "COPY INITIAL TEMPERATURE TO ELEMENTS NODES").buildAndExecute(*this);
////	CopyNodesToBoundaryNodes(*ParametersTemperature::Initial::output, temp.initial.boundary.node, "COPY INITIAL TEMPERATURE TO BOUNDARY NODES").buildAndExecute(*this);
////	ElementsGaussPointsBuilder<1>(integration.N, temp.initial.node, temp.initial.gp, "INTEGRATE INITIAL TEMPERATURE INTO ELEMENTS GAUSS POINTS").buildAndExecute(*this);
////	BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.initial.boundary.node, temp.initial.boundary.gp, "INTEGRATE INITIAL TEMPERATURE INTO BOUNDARY GAUSS POINTS").buildAndExecute(*this);
//
//	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end() || ParametersGradient::output) {
//		copyNodesToEnodes(*this, *ParametersTemperature::output, temp.node);
//	}
//
//	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end()) {
//		moveEnodesToGPs(*this, temp.node, temp.gp, 1);
//		Variable::list.egps["TEMPERATURE"] = new ParameterVariable(temp.gp.data, temp.gp.isconst, temp.gp.update, 0, 1);
//	}



//	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
//	CopyElementParameters(temp.initial.node, temp.node, "COPY INITIAL TEMPERATURE TO ELEMENT NODES").buildAndExecute(*this);
//	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node, "COPY TEMPERATURE TO ELEMENTS NODES"));
//	builders.push_back(new ElementsGaussPointsBuilder<1>(integration.N, temp.node, temp.gp, "INTEGRATE TEMPERATURE INTO ELEMENTS GAUSS POINTS"));
//	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node, "COPY TEMPERATURE TO BOUNDARY NODES"));
//	builders.push_back(new BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp, "INTEGRATE TEMPERATURE INTO BOUNDARY GAUSS POINTS"));

//	results();
	return correct;
}

void StructuralMechanics::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");
	initNames();

	bool correct = true;

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	_evaluate(); // fill coordinates, compute determinants
//	printElementVolume(integration.weight, integration.jacobiDeterminant);
//	printBoundarySurface(integration.boundary.weight, integration.boundary.jacobian);
	eslog::info(" ============================================================================================= \n");

	correct &= initTemperature();

	if (step::step.loadstep == 0) {
		///////////////////////////////////// Set materials and check if there is not any incorrect region intersection
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		eslog::info("\n  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");
		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			MaterialConfiguration *mat = info::mesh->materials[i];

			switch (mat->material_model) {
			case MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
				eslog::info("                                                                               LINEAR ELASTIC \n");
				if (info::mesh->dimension == 2) {
					switch (settings.element_behaviour) {
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
						eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRAIN \n");
						break;
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
						eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRESS \n");
						break;
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
						eslog::info("     ELEMENT BEHAVIOR:                                            PLANE STRESS WITH THICKNESS \n");
						correct &= examineElementParameter("THICKNESS", settings.thickness, thickness.gp.externalValues);
						fromExpression(*this, thickness.gp, thickness.gp.externalValues);
						break;
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
						eslog::info("     ELEMENT BEHAVIOR:                                                           AXISYMMETRIC \n");
						correct &= examineElementParameter("THICKNESS", settings.thickness, thickness.gp.externalValues);
						fromExpression(*this, thickness.gp, thickness.gp.externalValues);
						break;
					}
				}
				eslog::info("                                                                                               \n");

				switch (mat->linear_elastic_properties.model) {
				case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
					eslog::info("                MODEL:                                                              ISOTROPIC \n");
					correct &= examineMaterialParameter(mat->name, "EX", mat->linear_elastic_properties.young_modulus.get(0, 0), material.model.isoYoungModulus.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "mi", mat->linear_elastic_properties.poisson_ratio.get(0, 0), material.model.isoPoissonRatio.externalValues, 0);
					break;
				case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
					eslog::info("                MODEL:                                                            ORTHOTROPIC \n");
					break;
				case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
					eslog::info("                MODEL:                                                            ANISOTROPIC \n");
					break;
				}
				break;
			case MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
				eslog::info("                                                                                HYPER ELASTIC \n");
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValues, 0);
		}

//		fromExpression(*this, cooSystem.cartesian2D, cooSystem.cartesian2D.externalValues);
//		fromExpression(*this, cooSystem.cartesian3D, cooSystem.cartesian3D.externalValues);

		fromExpression(*this, material.model.isoYoungModulus, material.model.isoYoungModulus.externalValues);
		fromExpression(*this, material.model.isoPoissonRatio, material.model.isoPoissonRatio.externalValues);
		fromExpression(*this, material.model.youngModulus, material.model.youngModulus.externalValues);
		fromExpression(*this, material.model.poissonRatio, material.model.poissonRatio.externalValues);
		fromExpression(*this, material.model.shearModulus, material.model.shearModulus.externalValues);
		fromExpression(*this, material.model.anisotropic3D, material.model.anisotropic3D.externalValues);

		fromExpression(*this, material.density, material.density.externalValues);

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);

		elasticity(*this);

		eslog::info(" ============================================================================================= \n");
	}

	stiffness(*this);

	if (configuration.acceleration.size()) {
		correct &= examineElementParameter("ACCELERATION.X", configuration.acceleration, acceleration.gp.externalValues, 0);
		correct &= examineElementParameter("ACCELERATION.Y", configuration.acceleration, acceleration.gp.externalValues, 1);
		if (info::mesh->dimension == 3) {
			correct &= examineElementParameter("ACCELERATION.Z", configuration.acceleration, acceleration.gp.externalValues, 2);
		}
		fromExpression(*this, acceleration.gp, acceleration.gp.externalValues);
	}

	if (configuration.angular_velocity.size()) {
		switch (info::mesh->dimension) {
		case 3:
			correct &= examineElementParameter("ANGULAR_VELOCITY.X", configuration.angular_velocity, angularVevocity.gp.externalValues, 0);
			correct &= examineElementParameter("ANGULAR_VELOCITY.Y", configuration.angular_velocity, angularVevocity.gp.externalValues, 1);
			correct &= examineElementParameter("ANGULAR_VELOCITY.Z", configuration.angular_velocity, angularVevocity.gp.externalValues, 2);
			break;
		case 2:
			switch (settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				correct &= examineElementParameter("ANGULAR_VELOCITY.Z", configuration.angular_velocity, angularVevocity.gp.externalValues, 2);
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				correct &= examineElementParameter("ANGULAR_VELOCITY.Y", configuration.angular_velocity, angularVevocity.gp.externalValues, 1);
				break;
			}
		}

		fromExpression(*this, angularVevocity.gp, angularVevocity.gp.externalValues);
	}

	if (configuration.displacement.size()) {
		correct &= examineBoundaryParameter("DISPLACEMENT.X", configuration.displacement, displacement.node.externalValues, 0);
		correct &= examineBoundaryParameter("DISPLACEMENT.Y", configuration.displacement, displacement.node.externalValues, 1);
		if (info::mesh->dimension == 3) {
			correct &= examineBoundaryParameter("DISPLACEMENT.Z", configuration.displacement, displacement.node.externalValues, 2);
		}
		fromExpression(*this, displacement.node, displacement.node.externalValues);
	}

	if (configuration.normal_pressure.size()) {
		correct &= examineBoundaryParameter("NORMAL PRESSURE", configuration.normal_pressure, normalPressure.gp.externalValues);
		fromExpression(*this, normalPressure.gp, normalPressure.gp.externalValues);
	}
	RHS(*this);

	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void StructuralMechanics::connect(SteadyState &scheme)
{
	addFiller(*this, scheme);
}

void StructuralMechanics::evaluate(SteadyState &scheme)
{
	controller.setUpdate();
	reset(scheme.K, scheme.f, scheme.dirichlet);
//	printVersions();
	iterate();
	fill();
	update(scheme.K, scheme.f);
	controller.resetUpdate();
}

void StructuralMechanics::_evaluate()
{
	controller.setUpdate();
	iterate();
	controller.resetUpdate();
}

void StructuralMechanics::updateSolution(SteadyState &scheme)
{
	scheme.x->storeTo(Results::displacement->data);
	results(); // do we need an update mechanism?
}

void StructuralMechanics::initNames()
{
	integration.weight.name = "integration.weight";
	integration.N.name = "integration.N";
	integration.dN.name = "integration.dN";
	integration.dND.name = "integration.dND";
	integration.jacobiDeterminant.name = "integration.jacobiDeterminant";
	integration.jacobiInversion.name = "integration.jacobiInversion";

	coords.node.name = "coords.node";
	coords.gp.name = "coords.gp";

	thickness.gp.name = "thickness.gp";

	cooSystem.cartesian2D.name = "cooSystem.cartesian2D";
	cooSystem.cartesian3D.name = "cooSystem.cartesian3D";

	material.model.isoPoissonRatio.name = "material.model.isoPoissonRatio";
	material.model.isoYoungModulus.name = "material.model.isoYoungModulus";
	material.model.poissonRatio.name = "material.model.poissonRatio";
	material.model.youngModulus.name = "material.model.youngModulus";
	material.model.shearModulus.name = "material.model.shearModulus";
	material.model.anisotropic3D.name = "material.model.anisotropic3D";

	material.density.name = "material.density";
	material.mass.name = "material.mass";

	material.elasticityPlane.name = "material.elasticityPlane";
	material.elasticityAxisymm.name = "material.elasticityAxisymm";
	material.elasticity3D.name = "material.elasticity3D";

//	temp.initial.node.name = "temp.initial.node";
//	temp.initial.gp.name = "temp.initial.gp";
//	temp.node.name = "temp.node";
//	temp.gp.name = "temp.gp";
//
//	translationMotions.gp.name = "translationMotions.gp";
//	translationMotions.stiffness.name = "translationMotions.stiffness";
//	translationMotions.rhs.name = "translationMotions.rhs";
//
//	heatSource.gp.name = "heatSource.gp";

	elements.stiffness.name = "elements.stiffness";
//	elements.mass.name = "elements.mass";
	elements.rhs.name = "elements.rhs";

	acceleration.gp.name = "acceleration";
	angularVevocity.gp.name = "angularVevocity";

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		integration.boundary.N.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.N";
		integration.boundary.dN.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.dN";
		integration.boundary.normal.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.normal";
		integration.boundary.weight.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.weight";
		integration.boundary.jacobian.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.jacobian";

		coords.boundary.node.regions[r].name = info::mesh->boundaryRegions[r]->name + "::coords.boundary.node";
		coords.boundary.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::coords.boundary.gp";

		normalPressure.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::normalPressure.gp";
		displacement.node.regions[r].name = info::mesh->boundaryRegions[r]->name + "::displacement.node";

		elements.boundary.stiffness.regions[r].name = info::mesh->boundaryRegions[r]->name + "::stiffness";
		elements.boundary.rhs.regions[r].name = info::mesh->boundaryRegions[r]->name + "::rhs";
	}
}

void StructuralMechanics::printVersions()
{
	printParameterStats(integration.weight);
	printParameterStats(integration.N);
	printParameterStats(integration.dN);
	printParameterStats(integration.dND);
	printParameterStats(integration.jacobiDeterminant);
	printParameterStats(integration.jacobiInversion);

	printParameterStats(coords.node);
	printParameterStats(coords.gp);

	printParameterStats(thickness.gp);

	printParameterStats(cooSystem.cartesian2D);
	printParameterStats(cooSystem.cartesian3D);

	printParameterStats(material.model.isoPoissonRatio);
	printParameterStats(material.model.isoYoungModulus);
	printParameterStats(material.model.poissonRatio);
	printParameterStats(material.model.youngModulus);
	printParameterStats(material.model.shearModulus);
	printParameterStats(material.model.anisotropic3D);

	printParameterStats(material.density);
	printParameterStats(material.mass);

	printParameterStats(material.elasticityPlane);
	printParameterStats(material.elasticityAxisymm);
	printParameterStats(material.elasticity3D);

	printParameterStats(elements.stiffness);
//	printParameterStats("elements.mass", elements.mass);
	printParameterStats(elements.rhs);

	printParameterStats(acceleration.gp);
	printParameterStats(angularVevocity.gp);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		printf("REGION: %s\n", info::mesh->boundaryRegions[r]->name.c_str());

		printParameterStats(integration.boundary.N.regions[r]);
		printParameterStats(integration.boundary.dN.regions[r]);
		printParameterStats(integration.boundary.normal.regions[r]);
		printParameterStats(integration.boundary.weight.regions[r]);
		printParameterStats(integration.boundary.jacobian.regions[r]);

		printParameterStats(coords.boundary.node.regions[r]);
		printParameterStats(coords.boundary.gp.regions[r]);

		printParameterStats(normalPressure.gp.regions[r]);
		printParameterStats(displacement.node.regions[r]);

		printParameterStats(elements.boundary.stiffness.regions[r]);
		printParameterStats(elements.boundary.rhs.regions[r]);
	}
}



