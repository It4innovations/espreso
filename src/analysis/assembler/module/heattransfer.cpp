
#include "heattransfer.h"
#include "assembler.hpp"

#include "basis/expression/variable.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/assembler/operator.h"
#include "analysis/assembler/operators/operators.h"
#include "analysis/scheme/steadystate.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

NodeData* HeatTransfer::Results::temperature = nullptr;
NodeData* HeatTransfer::Results::initialTemperature = nullptr;
ElementData* HeatTransfer::Results::translationMotion = nullptr;
ElementData* HeatTransfer::Results::gradient = nullptr;
ElementData* HeatTransfer::Results::flux = nullptr;

HeatTransfer::HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{

}

void HeatTransfer::initParameters()
{
	if (Results::initialTemperature == nullptr) {
		Results::initialTemperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");

		Variable::list.node["INITIAL_TEMPERATURE"] = new OutputVariable(Results::initialTemperature, 0, 1);
		for (auto it = settings.initial_temperature.begin(); it != settings.initial_temperature.end(); ++it) {
			it->second.scope = ECFExpression::Scope::ENODES;
			for (auto p = it->second.parameters.begin(); p != it->second.parameters.end(); ++p) {
				Variable::list.enodes.insert(std::make_pair(*p, nullptr));
			}
		}
	}
	if (Results::temperature == nullptr) {
		Results::temperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
		Variable::list.node["TEMPERATURE"] = new OutputVariable(Results::temperature, 0, 1);
	}
	if (info::ecf->output.results_selection.translation_motions && Results::translationMotion == nullptr) {
		Results::translationMotion = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
	}
	if (info::ecf->output.results_selection.gradient && Results::gradient == nullptr) {
		Results::gradient = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
	}
	if (info::ecf->output.results_selection.flux && Results::flux == nullptr) {
		Results::flux = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
	}
}

bool HeatTransfer::initTemperature()
{
	// This code has to solve problem that initial temperature is set to elements regions, but we need it in nodes
	// 1. settings -> nodeInitialTemperature
	// 2. nodeInitialTemperature -> initialTemperature (here we average values)
	// 3. Dirichlet -> initialTemperature (if 'init_temp_respect_bc')
	// 4. initialTemperature -> nodeInitialTemperature (TODO: check the correction with TB)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool correct = true;

	if (configuration.temperature.size()) {
		correct &= examineBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature, temperature.node.externalValues);
		fromExpression(*this, temperature.node, temperature.node.externalValues);
	}

	if (settings.init_temp_respect_bc) {
		temp.initial.node.setConstness(false);
	}
	examineElementParameter("INITIAL TEMPERATURE", settings.initial_temperature, temp.initial.node.externalValues);
	// TODO: fix
//	fromExpression(*this, temp.initial.node, temp.initial.node.externalValues);
	_evaluate();

//	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
//		HeatTransfer::insertParameters(it->second.evaluator);
//	}
//
//	temp.initial.node.builder->buildAndExecute(*this);

	// TODO: fix
//	averageEnodesToNodes(temp.initial.node, *Results::initialTemperature);
//	Results::temperature->data = Results::initialTemperature->data;

//	if (info::mesh->dimension == 2 && info::ecf->heat_transfer_2d.init_temp_respect_bc) {
//		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::Initial::output, "SET INITIAL TEMPERATURE ACCORDING TO DIRICHLET").buildAndExecute(*this);
//	}
//	if (info::mesh->dimension == 3 && info::ecf->heat_transfer_3d.init_temp_respect_bc) {
//		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::Initial::output, "SET INITIAL TEMPERATURE ACCORDING TO DIRICHLET").buildAndExecute(*this);
//	}
//	CopyNodesToElementsNodes(*ParametersTemperature::Initial::output, temp.initial.node, "COPY INITIAL TEMPERATURE TO ELEMENTS NODES").buildAndExecute(*this);
//	CopyNodesToBoundaryNodes(*ParametersTemperature::Initial::output, temp.initial.boundary.node, "COPY INITIAL TEMPERATURE TO BOUNDARY NODES").buildAndExecute(*this);
//	ElementsGaussPointsBuilder<1>(integration.N, temp.initial.node, temp.initial.gp, "INTEGRATE INITIAL TEMPERATURE INTO ELEMENTS GAUSS POINTS").buildAndExecute(*this);
//	BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.initial.boundary.node, temp.initial.boundary.gp, "INTEGRATE INITIAL TEMPERATURE INTO BOUNDARY GAUSS POINTS").buildAndExecute(*this);

	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end() || Results::gradient) {
		copyNodesToEnodes(*this, *Results::temperature, temp.node);
	}

	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end()) {
		moveEnodesToGPs(*this, temp.node, temp.gp, 1);
		Variable::list.egps["TEMPERATURE"] = new ParameterVariable(temp.gp.data, temp.gp.isconst, temp.gp.update, 0, 1);
	}



//	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
//	CopyElementParameters(temp.initial.node, temp.node, "COPY INITIAL TEMPERATURE TO ELEMENT NODES").buildAndExecute(*this);
//	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node, "COPY TEMPERATURE TO ELEMENTS NODES"));
//	builders.push_back(new ElementsGaussPointsBuilder<1>(integration.N, temp.node, temp.gp, "INTEGRATE TEMPERATURE INTO ELEMENTS GAUSS POINTS"));
//	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node, "COPY TEMPERATURE TO BOUNDARY NODES"));
//	builders.push_back(new BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp, "INTEGRATE TEMPERATURE INTO BOUNDARY GAUSS POINTS"));

	results();
	return correct;
}

void HeatTransfer::analyze()
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
	printElementVolume(integration.weight, integration.jacobiDeterminant);
	printBoundarySurface(integration.boundary.weight, integration.boundary.jacobian);
	eslog::info(" ============================================================================================= \n");
	correct &= initTemperature();

//	Variable::print();

	if (step::step.loadstep == 0) {
		if (info::mesh->dimension == 2) {
			correct &= examineElementParameter("THICKNESS", settings.thickness, thickness.gp.externalValues);
			fromExpression2D(*this, thickness.gp, thickness.gp.externalValues, [] (auto &element) { return element.ecf.thickness; });
		}

		///////////////////////////////////// Set materials and check if there is not any incorrect region intersection
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		eslog::info("\n  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");
		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			MaterialConfiguration *mat = info::mesh->materials[i];

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				eslog::info("    COORDINATE SYSTEM:                                                              CARTESIAN \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian2D.externalValues, 0);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x, cooSystem.cartesian3D.externalValues, 0);
					examineMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y, cooSystem.cartesian3D.externalValues, 1);
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian3D.externalValues, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.spherical.externalValues, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.spherical.externalValues, 1);
					examineMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z, cooSystem.spherical.externalValues, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.externalValues, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.externalValues, 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.externalValues, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.externalValues, 1);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValues, 0);
			correct &= examineMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity, material.heatCapacity.externalValues, 0);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.isotropic.externalValues, 0);
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.externalValues, 1);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.externalValues, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.diagonal.externalValues, 2);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric2D.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric2D.externalValues, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric2D.externalValues, 1);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric3D.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric3D.externalValues, 3);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.symmetric3D.externalValues, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric3D.externalValues, 1);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.symmetric3D.externalValues, 4);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.symmetric3D.externalValues, 2);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.externalValues, 3);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.externalValues, 1);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.externalValues, 2);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.externalValues, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.externalValues, 4);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.anisotropic.externalValues, 8);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.externalValues, 1);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.anisotropic.externalValues, 5);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.anisotropic.externalValues, 2);
					correct &= examineMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.externalValues, 3);
					correct &= examineMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1), material.model.anisotropic.externalValues, 7);
					correct &= examineMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0), material.model.anisotropic.externalValues, 6);
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		fromExpression(*this, cooSystem.cartesian2D, cooSystem.cartesian2D.externalValues, [] (auto &element) { return element.ecf.angle; });
		fromExpression(*this, cooSystem.cartesian3D, cooSystem.cartesian3D.externalValues, [] (auto &element) { return element.ecf.angle; });
		fromExpression(*this, cooSystem.spherical, cooSystem.spherical.externalValues, [] (auto &element) { return element.ecf.angle; });
		fromExpression(*this, cooSystem.cylindric, cooSystem.cylindric.externalValues, [] (auto &element) { return element.ecf.angle; });

		fromExpression(*this, material.model.isotropic, material.model.isotropic.externalValues, [] (auto &element) { return element.conductivity; }); // it is possible to put it directly to target
		fromExpression(*this, material.model.diagonal, material.model.diagonal.externalValues, [] (auto &element) { return element.ecf.conductivity; });
		fromExpression(*this, material.model.symmetric2D, material.model.symmetric2D.externalValues, [] (auto &element) { return element.ecf.conductivity; });
		fromExpression(*this, material.model.symmetric3D, material.model.symmetric3D.externalValues, [] (auto &element) { return element.ecf.conductivity; });
		fromExpression(*this, material.model.anisotropic, material.model.anisotropic.externalValues, [] (auto &element) { return element.ecf.conductivity; });

		fromExpression(*this, material.density, material.density.externalValues, [] (auto &element) { return element.ecf.density; });
		fromExpression(*this, material.heatCapacity, material.heatCapacity.externalValues, [] (auto &element) { return element.ecf.heatCapacity; });

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);

		thermalConductivity(*this);

		eslog::info(" ============================================================================================= \n");
	}

	gradient.xi.resize(1);
	controller.prepare(gradient.xi);
	stiffness(*this);

	if (configuration.heat_source.size()) {
		correct &= examineElementParameter("HEAT SOURCE", configuration.heat_source, heatSource.gp.externalValues);
		fromExpression(*this, heatSource.gp, heatSource.gp.externalValues);
	}
	if (configuration.heat_flow.size()) {
		correct &= examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, heatFlow.gp.externalValues);
		fromExpression(*this, heatFlow.gp, heatFlow.gp.externalValues);
	}
	if (configuration.heat_flux.size()) {
		correct &= examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, heatFlux.gp.externalValues);
		fromExpression(*this, heatFlux.gp, heatFlux.gp.externalValues);
	}
	RHS(*this);

	outputGradient(*this);
	outputFlux(*this);

	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
	eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize<HeatTransfer, HeatTransferOperator>());
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void HeatTransfer::connect(SteadyState &scheme)
{
	addFiller(*this, scheme);
}

void HeatTransfer::evaluate(SteadyState &scheme)
{
	controller.setUpdate();
	reset(scheme.K, scheme.f, scheme.dirichlet);
//	iterate();
//	printVersions();
	eslog::info("       = SIMD LOOP                                                          %8.3f s = \n", assemble<HeatTransfer, HeatTransferOperator>());
//	fill();
	update(scheme.K, scheme.f);
	controller.resetUpdate();
}

void HeatTransfer::_evaluate()
{
	controller.setUpdate();
//	iterate();
//	printVersions();
	assemble<HeatTransfer, HeatTransferOperator>();
	controller.resetUpdate();
}

void HeatTransfer::updateSolution(SteadyState &scheme)
{
	scheme.x->storeTo(Results::temperature->data);
	results(); // do we need an update mechanism?
	temp.node.setUpdate(1);
}

void HeatTransfer::initNames()
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
	cooSystem.cylindric.name = "cooSystem.cylindric";
	cooSystem.spherical.name = "cooSystem.spherical";

	material.model.isotropic.name = "material.model.isotropic";
	material.model.diagonal.name = "material.model.diagonal";
	material.model.symmetric2D.name = "material.model.symmetric2D";
	material.model.symmetric3D.name = "material.model.symmetric3D";
	material.model.anisotropic.name = "material.model.anisotropic";

	material.conductivityIsotropic.name = "material.conductivityIsotropic";
	material.conductivity.name = "material.conductivity";
	material.density.name = "material.density";
	material.heatCapacity.name = "material.heatCapacity";
	material.mass.name = "material.mass";

	temp.initial.node.name = "temp.initial.node";
	temp.initial.gp.name = "temp.initial.gp";
	temp.node.name = "temp.node";
	temp.gp.name = "temp.gp";

	translationMotions.gp.name = "translationMotions.gp";
	translationMotions.stiffness.name = "translationMotions.stiffness";
	translationMotions.rhs.name = "translationMotions.rhs";

	heatSource.gp.name = "heatSource.gp";

	elements.stiffness.name = "elements.stiffness";
	elements.mass.name = "elements.mass";
	elements.rhs.name = "elements.rhs";

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		integration.boundary.N.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.N";
		integration.boundary.dN.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.dN";
		integration.boundary.weight.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.weight";
		integration.boundary.jacobian.regions[r].name = info::mesh->boundaryRegions[r]->name + "::integration.boundary.jacobian";

		coords.boundary.node.regions[r].name = info::mesh->boundaryRegions[r]->name + "::coords.boundary.node";
		coords.boundary.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::coords.boundary.gp";

		convection.heatTransferCoeficient.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::convection.heatTransferCoeficient.gp";
		convection.externalTemperature.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::convection.externalTemperature.gp";

		heatFlow.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::heatFlow.gp";
		heatFlux.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::heatFlux.gp";
		q.gp.regions[r].name = info::mesh->boundaryRegions[r]->name + "::q.gp";
	}
}

void HeatTransfer::printVersions()
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
	printParameterStats(cooSystem.cylindric);
	printParameterStats(cooSystem.spherical);

	printParameterStats(material.model.isotropic);
	printParameterStats(material.model.diagonal);
	printParameterStats(material.model.symmetric2D);
	printParameterStats(material.model.symmetric3D);
	printParameterStats(material.model.anisotropic);

	printParameterStats(material.conductivityIsotropic);
	printParameterStats(material.conductivity);
	printParameterStats(material.density);
	printParameterStats(material.heatCapacity);
	printParameterStats(material.mass);

	printParameterStats(temp.initial.node);
	printParameterStats(temp.initial.gp);
	printParameterStats(temp.node);
	printParameterStats(temp.gp);

	printParameterStats(translationMotions.gp);
	printParameterStats(translationMotions.stiffness);
	printParameterStats(translationMotions.rhs);

	printParameterStats(heatSource.gp);

	printParameterStats(elements.stiffness);
	printParameterStats(elements.mass);
	printParameterStats(elements.rhs);

	printParameterStats(gradient.xi);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		printf("REGION: %s\n", info::mesh->boundaryRegions[r]->name.c_str());

		printParameterStats(integration.boundary.N.regions[r]);
		printParameterStats(integration.boundary.dN.regions[r]);
		printParameterStats(integration.boundary.weight.regions[r]);
		printParameterStats(integration.boundary.jacobian.regions[r]);

		printParameterStats(coords.boundary.node.regions[r]);
		printParameterStats(coords.boundary.gp.regions[r]);

		printParameterStats(convection.heatTransferCoeficient.gp.regions[r]);
		printParameterStats(convection.externalTemperature.gp.regions[r]);

		printParameterStats(heatFlow.gp.regions[r]);
		printParameterStats(heatFlux.gp.regions[r]);
		printParameterStats(q.gp.regions[r]);
	}
}
