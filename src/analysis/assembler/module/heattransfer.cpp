
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
: gsettings(gsettings), configuration(configuration), K{}, M{}, rhs{}, x{}
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

		Variable::list.node["INITIAL_TEMPERATURE"] = Variable(0, 1, ParametersTemperature::Initial::output->data.data(), false, true);
		for (auto it = initialTemperature.begin(); it != initialTemperature.end(); ++it) {
			it->second.scope = ECFExpression::Scope::ENODES;
			for (auto p = it->second.parameters.begin(); p != it->second.parameters.end(); ++p) {
				Variable::list.enodes.insert(std::make_pair(*p, Variable()));
			}
		}
	}
	if (ParametersTemperature::output == nullptr) {
		ParametersTemperature::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
		Variable::list.node["TEMPERATURE"] = Variable(0, 1, ParametersTemperature::output->data.data(), false, true);
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

void AX_HeatTransfer::initTemperature()
{
	// This code has to solve problem that initial temperature is set to elements regions, but we need it in nodes
	// 1. settings -> nodeInitialTemperature
	// 2. nodeInitialTemperature -> initialTemperature (here we average values)
	// 3. Dirichlet -> initialTemperature (if 'init_temp_respect_bc')
	// 4. initialTemperature -> nodeInitialTemperature (TODO: check the correction with TB)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
//	iterate();
	evaluateFromExpression(*this, temp.initial.node, temp.initial.node.externalValue);

//	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
//		AX_HeatTransfer::insertParameters(it->second.evaluator);
//	}
//
//	temp.initial.node.builder->buildAndExecute(*this);

	averageEnodesToNodes(temp.initial.node, *ParametersTemperature::Initial::output);

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

	auto it = Variable::list.egps.end();
	if ((it = Variable::list.egps.find("TEMPERATURE"))!= Variable::list.egps.end()) {
		copyNodesToEnodes(*this, *ParametersTemperature::output, temp.node);
		moveEnodesToGPs(*this, temp.node, temp.gp, 1);
		Variable::list.egps["TEMPERATURE"] = Variable(0, 1, temp.gp.data->datatarray().data(), false, true);
	}
//	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
//	CopyElementParameters(temp.initial.node, temp.node, "COPY INITIAL TEMPERATURE TO ELEMENT NODES").buildAndExecute(*this);
//	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node, "COPY TEMPERATURE TO ELEMENTS NODES"));
//	builders.push_back(new ElementsGaussPointsBuilder<1>(integration.N, temp.node, temp.gp, "INTEGRATE TEMPERATURE INTO ELEMENTS GAUSS POINTS"));
//	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node, "COPY TEMPERATURE TO BOUNDARY NODES"));
//	builders.push_back(new BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp, "INTEGRATE TEMPERATURE INTO BOUNDARY GAUSS POINTS"));
}


void AX_HeatTransfer::init(AX_SteadyState &scheme)
{
	K = scheme.K;
	rhs = scheme.f;
	x = scheme.x;

	analyze();
}

void AX_HeatTransfer::analyze()
{
	eslog::info(" == PHYSICS                                                                   HEAT TRANSFER == \n");
	eslog::info(" ============================================================================================= \n");
	bool correct = true;

	if (info::mesh->dimension == 2) {
		setMaterials(info::ecf->heat_transfer_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_2d.material_set);
		validateRegionSettings("THICKNESS", info::ecf->heat_transfer_2d.thickness);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature);
	}
	if (info::mesh->dimension == 3) {
		setMaterials(info::ecf->heat_transfer_3d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_3d.material_set);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature);
	}

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);

	initTemperature();

	elementIntegration(*this);

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
				eslog::info("         CONDUCTIVITY:                                                             SYMMETRIC \n");
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
	heatStiffness(*this);

	if (configuration.temperature.size()) {
		correct &= examineBoundaryParameter("TEMPERATURE", configuration.temperature, dirichlet.gp.externalValues);
		fromExpression(*this, dirichlet.gp, dirichlet.gp.externalValues);
	}
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

	eslog::info(" ============================================================================================= \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	} else {
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       FAIL  \n");
	}
	eslog::info(" ============================================================================================= \n\n");
	if (!correct) {
		eslog::globalerror("                                                               INVALID CONFIGURATION DETECTED \n");
	}
}

void AX_HeatTransfer::next()
{
	updateVersions();

	if (K != nullptr) {
		K->fill(0);
		K->touched = true;
	}
	if (M != nullptr) {
		M->fill(0);
		M->touched = true;
	}
	if (rhs != nullptr) {
		rhs->fill(0);
		rhs->touched = true;
	}

	iterate();

	printVersions();
}

void AX_HeatTransfer::initDirichlet(Vector_Sparse<double> &dirichlet)
{
	size_t dsize = 0;
	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dsize += region->nodes->datatarray().size();
	}
	dirichletIndices.reserve(dsize);
	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dirichletIndices.insert(dirichletIndices.end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
	}
	dirichletPermutation.resize(dsize);
	std::iota(dirichletPermutation.begin(), dirichletPermutation.end(), 0);
	std::sort(dirichletPermutation.begin(), dirichletPermutation.end(), [&] (const esint &i, const esint &j) { return dirichletIndices[i] < dirichletIndices[j]; });
	dsize = 0;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			++dsize;
		}
	}
	dirichlet.resize(info::mesh->nodes->IDs->datatarray().size(), dsize);
	auto dir = dirichlet.indices;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			*dir++ = dirichletIndices[*i];
		}
	}
}

void AX_HeatTransfer::fillDirichlet(Vector_Sparse<double> &dirichlet)
{
	size_t offset = 0;
	std::vector<double> values(dirichletPermutation.size());
	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelectedSparse(
				region->nodes->datatarray().size(),
				region->nodes->datatarray().data(),
				it->second.evaluator->params,
				values.data() + offset);
		offset += region->nodes->datatarray().size();
	}

	for (size_t i = 0; i < dirichletPermutation.size(); ++i) {
		dirichlet.vals[i] = values[dirichletPermutation[i]];
	}
	dirichlet.touched = true;
}

void AX_HeatTransfer::updateSolution()
{
	x->store(ParametersTemperature::output->data);
}

void AX_HeatTransfer::printVersions()
{
	printParamterStats("integration.weight", integration.weight);
	printParamterStats("integration.N", integration.N);
	printParamterStats("integration.dN", integration.dN);
	printParamterStats("integration.dND", integration.dND);
	printParamterStats("integration.jacobiDeterminant", integration.jacobiDeterminant);
	printParamterStats("integration.jacobiInversion", integration.jacobiInversion);

	printParamterStats("coords.node", coords.node);
	printParamterStats("coords.gp", coords.gp);

	printParamterStats("thickness.gp", thickness.gp);

	printParamterStats("cooSystem.cartesian2D", cooSystem.cartesian2D);
	printParamterStats("cooSystem.cartesian3D", cooSystem.cartesian3D);
	printParamterStats("cooSystem.cylindric", cooSystem.cylindric);
	printParamterStats("cooSystem.spherical", cooSystem.spherical);

	printParamterStats("material.model.isotropic", material.model.isotropic);
	printParamterStats("material.model.diagonal", material.model.diagonal);
	printParamterStats("material.model.symmetric2D", material.model.symmetric2D);
	printParamterStats("material.model.symmetric3D", material.model.symmetric3D);
	printParamterStats("material.model.anisotropic", material.model.anisotropic);

	printParamterStats("material.conductivityIsotropic", material.conductivityIsotropic);
	printParamterStats("material.conductivity", material.conductivity);
	printParamterStats("material.density", material.density);
	printParamterStats("material.heatCapacity", material.heatCapacity);
	printParamterStats("material.mass", material.mass);

	printParamterStats("temp.initial.output", temp.initial.output);
	printParamterStats("temp.initial.node", temp.initial.node);
	printParamterStats("temp.initial.gp", temp.initial.gp);
	printParamterStats("temp.output", temp.output);
	printParamterStats("temp.node", temp.node);
	printParamterStats("temp.gp", temp.gp);

	printParamterStats("translationMotions.output", translationMotions.output);
	printParamterStats("translationMotions.gp", translationMotions.gp);
	printParamterStats("translationMotions.stiffness", translationMotions.stiffness);
	printParamterStats("translationMotions.rhs", translationMotions.rhs);

	printParamterStats("elements.stiffness", elements.stiffness);
	printParamterStats("elements.mass", elements.mass);
	printParamterStats("elements.rhs", elements.rhs);

	if (gradient.output)
	{
		printParamterStats("gradient.output", gradient.output);
	}

	printParamterStats("gradient.xi", gradient.xi);

	if (flux.output)
	{
		printParamterStats("flux.output", flux.output);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		printf("REGION: %s\n", info::mesh->boundaryRegions[r]->name.c_str());
		printParamterStats("convection.heatTransferCoeficient.gp", convection.heatTransferCoeficient.gp.regions[r]);
		printParamterStats("convection.externalTemperature.gp", convection.externalTemperature.gp.regions[r]);

		printParamterStats("heatFlow.gp", heatFlow.gp.regions[r]);
		printParamterStats("heatFlux.gp", heatFlux.gp.regions[r]);
		printParamterStats("q.gp", q.gp.regions[r]);
	}
}
