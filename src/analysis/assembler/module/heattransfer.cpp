
#include "heattransfer.h"
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

AX_HeatTransfer::AX_HeatTransfer(AX_HeatTransfer *previous, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), K{}, M{}, rhs{}, x{}, dirichlet{}
{

}

void AX_HeatTransfer::initParameters()
{
	if (info::mesh->dimension == 2) {
		initialTemperature = info::ecf->heat_transfer_2d.initial_temperature;
	}
	if (info::mesh->dimension == 3) {
		initialTemperature = info::ecf->heat_transfer_3d.initial_temperature;
	}

	if (ParametersTemperature::Initial::output == nullptr) {
		ParametersTemperature::Initial::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");

		Variable::list.node["INITIAL_TEMPERATURE"] = new OutputVariable(ParametersTemperature::Initial::output, 0, 1);
		for (auto it = initialTemperature.begin(); it != initialTemperature.end(); ++it) {
			it->second.scope = ECFExpression::Scope::ENODES;
			for (auto p = it->second.parameters.begin(); p != it->second.parameters.end(); ++p) {
				Variable::list.enodes.insert(std::make_pair(*p, nullptr));
			}
		}
	}
	if (ParametersTemperature::output == nullptr) {
		ParametersTemperature::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
		Variable::list.node["TEMPERATURE"] = new OutputVariable(ParametersTemperature::output, 0, 1);
	}
	if (info::ecf->output.results_selection.translation_motions && ParametersTranslationMotions::output == nullptr) {
		ParametersTranslationMotions::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
	}
	if (info::ecf->output.results_selection.gradient && ParametersGradient::output == nullptr) {
		ParametersGradient::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
	}
	if (info::ecf->output.results_selection.flux && ParametersFlux::output == nullptr) {
		ParametersFlux::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
	}
}

bool AX_HeatTransfer::initTemperature()
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

	if (info::mesh->dimension == 2) {
		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
			temp.initial.node.setConstness(false);
		}
		examineElementParameter("INITIAL TEMPERATURE", initialTemperature, temp.initial.node.externalValue);
	}
	if (info::mesh->dimension == 3) {
		if (info::ecf->heat_transfer_3d.init_temp_respect_bc) {
			temp.initial.node.setConstness(false);
		}
		examineElementParameter("INITIAL TEMPERATURE", initialTemperature, temp.initial.node.externalValue);
	}
	fromExpression(*this, temp.initial.node, temp.initial.node.externalValue);
	_evaluate();

//	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
//		AX_HeatTransfer::insertParameters(it->second.evaluator);
//	}
//
//	temp.initial.node.builder->buildAndExecute(*this);

	averageEnodesToNodes(temp.initial.node, *ParametersTemperature::Initial::output);
	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;

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

	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end() || ParametersGradient::output) {
		copyNodesToEnodes(*this, *ParametersTemperature::output, temp.node);
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

void AX_HeatTransfer::printVolume()
{
	std::vector<double> volume(info::mesh->elementsRegions.size());
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			if (info::mesh->elements->eintervals[i].region == (esint)r || (info::mesh->elements->eintervals[i].region == 0 && r == info::mesh->elementsRegions.size() - 1)) {
				auto begin = (integration.jacobiDeterminant.data->begin() + i)->data();
				auto end = (integration.jacobiDeterminant.data->begin() + (i + 1))->data();
				for (auto det = begin; det != end; ++det) {
					volume[0] += *det;
					volume[r] += *det;
				}
			}
		}
	}

	std::vector<double> surface(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			auto begin = integration.boundary.jacobian.regions[r].data->datatarray().begin();
			auto end = integration.boundary.jacobian.regions[r].data->datatarray().end();
			for (auto det = begin; det != end; ++det) {
				surface[r] += *det;
			}
		}
	}

	eslog::info("  ELEMENT REGION VOLUME                                                                        \n");
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
	for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
		eslog::info("     %30s :                                            %e   \n", info::mesh->elementsRegions[r]->name.c_str(), volume[r]);
	}
	eslog::info("\n  BOUDNARY REGION SURFACE                                                                      \n");
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			eslog::info("     %30s :                                            %e   \n", info::mesh->boundaryRegions[r]->name.c_str(), surface[r]);
			info::mesh->boundaryRegions[r]->area = surface[r];
		}
	}
	eslog::info(" ============================================================================================= \n");
}

void AX_HeatTransfer::init(AX_SteadyState &scheme)
{
	this->K = scheme.K;
	this->rhs = scheme.f;
	this->x = scheme.x;
	this->dirichlet = scheme.dirichlet;

	initNames();
	analyze();
}

void AX_HeatTransfer::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");
	bool correct = true;

	if (info::mesh->dimension == 2) {
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_2d.material_set);
		validateRegionSettings("THICKNESS", info::ecf->heat_transfer_2d.thickness);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature);
	}
	if (info::mesh->dimension == 3) {
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_3d.material_set);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature);
	}

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	_evaluate(); // fill coordinates, compute determinants
	printVolume();

	correct &= initTemperature();

//	Variable::print();

	if (step::step.loadstep == 0) {
		if (info::mesh->dimension == 2) {
			correct &= examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, thickness.gp.externalValue);
			fromExpression(*this, thickness.gp, thickness.gp.externalValue);
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
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian2D.externalValue, 0);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x, cooSystem.cartesian3D.externalValue, 0);
					examineMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y, cooSystem.cartesian3D.externalValue, 1);
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian3D.externalValue, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.spherical.externalValue, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.spherical.externalValue, 1);
					examineMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z, cooSystem.spherical.externalValue, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.externalValue, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.externalValue, 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.externalValue, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.externalValue, 1);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValue, 0);
			correct &= examineMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity, material.heatCapacity.externalValue, 0);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.isotropic.externalValue, 0);
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.externalValue, 1);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.externalValue, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.diagonal.externalValue, 2);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric2D.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric2D.externalValue, 1);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric2D.externalValue, 2);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric3D.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric3D.externalValue, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.symmetric3D.externalValue, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric3D.externalValue, 3);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.symmetric3D.externalValue, 4);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.symmetric3D.externalValue, 5);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.externalValue, 1);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.externalValue, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.externalValue, 3);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.externalValue, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.externalValue, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.anisotropic.externalValue, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.externalValue, 3);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.anisotropic.externalValue, 4);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.anisotropic.externalValue, 5);
					correct &= examineMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.externalValue, 6);
					correct &= examineMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1), material.model.anisotropic.externalValue, 7);
					correct &= examineMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0), material.model.anisotropic.externalValue, 8);
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		fromExpression(*this, cooSystem.cartesian2D, cooSystem.cartesian2D.externalValue);
		fromExpression(*this, cooSystem.cartesian3D, cooSystem.cartesian3D.externalValue);
		fromExpression(*this, cooSystem.spherical, cooSystem.spherical.externalValue);
		fromExpression(*this, cooSystem.cylindric, cooSystem.cylindric.externalValue);

		fromExpression(*this, material.model.isotropic, material.model.isotropic.externalValue);
		fromExpression(*this, material.model.diagonal, material.model.diagonal.externalValue);
		fromExpression(*this, material.model.symmetric2D, material.model.symmetric2D.externalValue);
		fromExpression(*this, material.model.symmetric3D, material.model.symmetric3D.externalValue);
		fromExpression(*this, material.model.anisotropic, material.model.anisotropic.externalValue);

		fromExpression(*this, material.density, material.density.externalValue);
		fromExpression(*this, material.heatCapacity, material.heatCapacity.externalValue);

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->heat_transfer_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->heat_transfer_3d.material_set);
		}

		thermalConductivity(*this);

		eslog::info(" ============================================================================================= \n");
	}

	gradient.xi.resize(1);
	controller.prepare(gradient.xi);
	heatStiffness(*this);

	if (configuration.heat_flow.size()) {
		correct &= examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, heatFlow.gp.externalValues);
		fromExpression(*this, heatFlow.gp, heatFlow.gp.externalValues);
	}
	if (configuration.heat_flux.size()) {
		correct &= examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, heatFlux.gp.externalValues);
		fromExpression(*this, heatFlux.gp, heatFlux.gp.externalValues);
	}
	heatRHS(*this);

	addFiller(*this);

	outputGradient(*this);
	outputFlux(*this);

	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void AX_HeatTransfer::evaluate()
{
	controller.setUpdate();
	reset(K, M, rhs, dirichlet);
	iterate();
	fill();
	update(K, M, rhs);
	controller.resetUpdate();
}

void AX_HeatTransfer::_evaluate()
{
	controller.setUpdate();
	iterate();
	controller.resetUpdate();
}

void AX_HeatTransfer::updateSolution()
{
	x->store(ParametersTemperature::output->data);
	results(); // do we need an update mechanism?
	temp.node.setUpdate(1);
}

void AX_HeatTransfer::initNames()
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

void AX_HeatTransfer::printVersions()
{
	printParameterStats("integration.weight", integration.weight);
	printParameterStats("integration.N", integration.N);
	printParameterStats("integration.dN", integration.dN);
	printParameterStats("integration.dND", integration.dND);
	printParameterStats("integration.jacobiDeterminant", integration.jacobiDeterminant);
	printParameterStats("integration.jacobiInversion", integration.jacobiInversion);

	printParameterStats("coords.node", coords.node);
	printParameterStats("coords.gp", coords.gp);

	printParameterStats("thickness.gp", thickness.gp);

	printParameterStats("cooSystem.cartesian2D", cooSystem.cartesian2D);
	printParameterStats("cooSystem.cartesian3D", cooSystem.cartesian3D);
	printParameterStats("cooSystem.cylindric", cooSystem.cylindric);
	printParameterStats("cooSystem.spherical", cooSystem.spherical);

	printParameterStats("material.model.isotropic", material.model.isotropic);
	printParameterStats("material.model.diagonal", material.model.diagonal);
	printParameterStats("material.model.symmetric2D", material.model.symmetric2D);
	printParameterStats("material.model.symmetric3D", material.model.symmetric3D);
	printParameterStats("material.model.anisotropic", material.model.anisotropic);

	printParameterStats("material.conductivityIsotropic", material.conductivityIsotropic);
	printParameterStats("material.conductivity", material.conductivity);
	printParameterStats("material.density", material.density);
	printParameterStats("material.heatCapacity", material.heatCapacity);
	printParameterStats("material.mass", material.mass);

	printParameterStats("temp.initial.output", temp.initial.output);
	printParameterStats("temp.initial.node", temp.initial.node);
	printParameterStats("temp.initial.gp", temp.initial.gp);
	printParameterStats("temp.output", temp.output);
	printParameterStats("temp.node", temp.node);
	printParameterStats("temp.gp", temp.gp);

	printParameterStats("translationMotions.output", translationMotions.output);
	printParameterStats("translationMotions.gp", translationMotions.gp);
	printParameterStats("translationMotions.stiffness", translationMotions.stiffness);
	printParameterStats("translationMotions.rhs", translationMotions.rhs);

	printParameterStats("elements.stiffness", elements.stiffness);
	printParameterStats("elements.mass", elements.mass);
	printParameterStats("elements.rhs", elements.rhs);

	if (gradient.output)
	{
		printParameterStats("gradient.output", gradient.output);
	}

	printParameterStats("gradient.xi", gradient.xi);

	if (flux.output)
	{
		printParameterStats("flux.output", flux.output);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		printf("REGION: %s\n", info::mesh->boundaryRegions[r]->name.c_str());

		printParameterStats("integration.boundary.N", integration.boundary.N.regions[r]);
		printParameterStats("integration.boundary.dN", integration.boundary.dN.regions[r]);
		printParameterStats("integration.boundary.weight", integration.boundary.weight.regions[r]);
		printParameterStats("integration.boundary.jacobian", integration.boundary.jacobian.regions[r]);

		printParameterStats("coords.boundary.node", coords.boundary.node.regions[r]);
		printParameterStats("coords.boundary.gp", coords.boundary.gp.regions[r]);

		printParameterStats("convection.heatTransferCoeficient.gp", convection.heatTransferCoeficient.gp.regions[r]);
		printParameterStats("convection.externalTemperature.gp", convection.externalTemperature.gp.regions[r]);

		printParameterStats("heatFlow.gp", heatFlow.gp.regions[r]);
		printParameterStats("heatFlux.gp", heatFlux.gp.regions[r]);
		printParameterStats("q.gp", q.gp.regions[r]);
	}
}
