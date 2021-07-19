
#include "heattransfer.h"
#include "assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"

#include "analysis/assembler/operators/operators.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

#include "basis/utilities/print.h"

using namespace espreso;

AX_HeatTransfer::AX_HeatTransfer(AX_HeatTransfer *previous, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration)
{

}

void AX_HeatTransfer::initTemperature()
{
	/// This code has to solve problem that initial temperature is set to elements regions, but we need it in nodes
	/// 1. settings -> nodeInitialTemperature
	/// 2. nodeInitialTemperature -> initialTemperature (here we average values)
	/// 3. Dirichlet -> initialTemperature (if 'init_temp_respect_bc')
	/// 4. initialTemperature -> nodeInitialTemperature (TODO: check the correction with TB)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	temp.initial.node.builder = new ExpressionsToElementsSimple(temp.initial.node, 273.15, "INITIAL TEMPERATURE");
//	if (info::mesh->dimension == 2) {
//		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
//			temp.initial.node.setConstness(false);
//		}
//		examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature, *temp.initial.node.builder);
//	}
//	if (info::mesh->dimension == 3) {
//		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
//			temp.initial.node.setConstness(false);
//		}
//		examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature, *temp.initial.node.builder);
//	}
//	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
////		AX_HeatTransfer::insertParameters(it->second.evaluator);
//	}
//
//	temp.initial.node.builder->buildAndExecute(*this);
//	AverageElementsNodesToNodes(temp.initial.node, *ParametersTemperature::Initial::output, "AVERAGE NODE INITIAL TEMPERATURE").buildAndExecute(*this);
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
//
//	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
//	CopyElementParameters(temp.initial.node, temp.node, "COPY INITIAL TEMPERATURE TO ELEMENT NODES").buildAndExecute(*this);
//	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node, "COPY TEMPERATURE TO ELEMENTS NODES"));
//	builders.push_back(new ElementsGaussPointsBuilder<1>(integration.N, temp.node, temp.gp, "INTEGRATE TEMPERATURE INTO ELEMENTS GAUSS POINTS"));
//	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node, "COPY TEMPERATURE TO BOUNDARY NODES"));
//	builders.push_back(new BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp, "INTEGRATE TEMPERATURE INTO BOUNDARY GAUSS POINTS"));
}


void AX_HeatTransfer::init()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info("  PHYSICS                                                                    HEAT TRANSFER 2D  \n");
	eslog::info(" ============================================================================================= \n");

	ParametersTemperature::Initial::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
	ParametersTemperature::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");

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

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	bool correct = true;

	if (step::step.loadstep == 0) {
		if (info::mesh->dimension == 2) {
			correct &= examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, thickness.gp.evaluator);
			fromExpression(*this, thickness.gp, thickness.gp.evaluator);
		}

		///////////////////////////////////// Set materials and check if there is not any incorrect region intersection
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		eslog::info("  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");
		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			const MaterialConfiguration *mat = info::mesh->materials[i];

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				eslog::info("    COORDINATE SYSTEM:                                                              CARTESIAN \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian2D.evaluator, 0);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x, cooSystem.cartesian3D.evaluator, 0);
					examineMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y, cooSystem.cartesian3D.evaluator, 1);
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, cooSystem.cartesian3D.evaluator, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.spherical.evaluator, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.spherical.evaluator, 1);
					examineMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z, cooSystem.spherical.evaluator, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.evaluator, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.evaluator, 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, cooSystem.cylindric.evaluator, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, cooSystem.cylindric.evaluator, 1);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.evaluator, 0);
			correct &= examineMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity, material.heatCapacity.evaluator, 0);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.isotropic.evaluator, 0);
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.evaluator, 1);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.diagonal.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.diagonal.evaluator, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.diagonal.evaluator, 2);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                             SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric2D.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric2D.evaluator, 1);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric2D.evaluator, 2);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.symmetric3D.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.symmetric3D.evaluator, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.symmetric3D.evaluator, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.symmetric3D.evaluator, 3);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.symmetric3D.evaluator, 4);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.symmetric3D.evaluator, 5);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.evaluator, 1);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.evaluator, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.evaluator, 3);
				}
				if (info::mesh->dimension == 3) {
					correct &= examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), material.model.anisotropic.evaluator, 0);
					correct &= examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), material.model.anisotropic.evaluator, 1);
					correct &= examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), material.model.anisotropic.evaluator, 2);
					correct &= examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), material.model.anisotropic.evaluator, 3);
					correct &= examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), material.model.anisotropic.evaluator, 4);
					correct &= examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), material.model.anisotropic.evaluator, 5);
					correct &= examineMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0), material.model.anisotropic.evaluator, 6);
					correct &= examineMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1), material.model.anisotropic.evaluator, 7);
					correct &= examineMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0), material.model.anisotropic.evaluator, 8);
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		fromExpression(*this, cooSystem.cartesian2D, cooSystem.cartesian2D.evaluator);
		fromExpression(*this, cooSystem.cartesian3D, cooSystem.cartesian3D.evaluator);
		fromExpression(*this, cooSystem.spherical, cooSystem.spherical.evaluator);
		fromExpression(*this, cooSystem.cylindric, cooSystem.cylindric.evaluator);

		fromExpression(*this, material.model.isotropic, material.model.isotropic.evaluator);
		fromExpression(*this, material.model.diagonal, material.model.diagonal.evaluator);
		fromExpression(*this, material.model.symmetric2D, material.model.symmetric2D.evaluator);
		fromExpression(*this, material.model.symmetric3D, material.model.symmetric3D.evaluator);
		fromExpression(*this, material.model.anisotropic, material.model.anisotropic.evaluator);

		fromExpression(*this, material.density, material.density.evaluator);
		fromExpression(*this, material.heatCapacity, material.heatCapacity.evaluator);

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

	eslog::info(" ============================================================================================= \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	} else {
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       FAIL  \n");
	}
	eslog::info(" ============================================================================================= \n");
	if (!correct) {
		eslog::globalerror("                                                               INVALID CONFIGURATION DETECTED \n");
	}
}

void AX_HeatTransfer::next()
{
	updateVersions();

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				size_t elementsInInterval = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;

				// First pass through operators
				for(size_t element = 0; element < elementsInInterval; ++element) {
					for (auto op = actionOps[i].begin(); op != actionOps[i].end(); ++op) {
						if((*op)->update) {
							if(element == 0 || !(*op)->isconst) {
								(**op)();
								++(**op);
							}
						}
					}
				}
			}
		}
	}

	std::cout << "COO[nd]: " << *coords.node.data << "\n";
	std::cout << "COO[gp]: " << *coords.gp.data << "\n";
	std::cout << "thick: " << *thickness.gp.data << "\n";
	std::cout << "cart: " << *cooSystem.cartesian2D.data << "\n";
	std::cout << "iso: " << *material.model.isotropic.data << "\n";
	std::cout << "dia: " << *material.model.diagonal.data << "\n";
	std::cout << "cond-iso: " << *material.conductivityIsotropic.data << "\n";
	std::cout << "cond-dia: " << *material.conductivity.data << "\n";
	std::cout << "stiffness: " << *elements.stiffness.data << "\n";
}

void AX_HeatTransfer::fillDirichletIndices(Vector_Sparse<double> &dirichlet)
{
	size_t dsize = 0;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *bregion = info::mesh->boundaryRegions[r];
		auto it = configuration.temperature.find(bregion->name);
		if (it != configuration.temperature.end()) {
			dsize += bregion->nodes->datatarray().size();
		}
	}
	dirichletIndices.reserve(dsize);
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *bregion = info::mesh->boundaryRegions[r];
		auto it = configuration.temperature.find(bregion->name);
		if (it != configuration.temperature.end()) {
			for (auto n = bregion->nodes->datatarray().cbegin(); n != bregion->nodes->datatarray().cend(); ++n) {
				dirichletIndices.push_back(*n);
			}
		}
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

bool AX_HeatTransfer::fillK(Matrix_Base<double> *K, const ElementMapping<double> &mapping)
{
	K->fill(0);
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
//		MatricesFiller(*this, mapping.elements[i].data, mapping.elements[i].position).apply(i);
	}
	return true;
}

bool AX_HeatTransfer::fillRHS(Vector_Base<double> *rhs, const ElementMapping<double> &mapping)
{
	rhs->fill(0);
	return true;
}

bool AX_HeatTransfer::fillDirichlet(Vector_Sparse<double> &dirichlet)
{
	std::vector<double> values;
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *bregion = info::mesh->boundaryRegions[r];
		auto it = configuration.temperature.find(bregion->name);
		if (it != configuration.temperature.end()) {
			for (auto n = bregion->nodes->datatarray().cbegin(); n != bregion->nodes->datatarray().cend(); ++n) {
				values.push_back(r);
			}
		}
	}
	for (size_t i = 0; i < dirichletPermutation.size(); ++i) {
		dirichlet.vals[i] = values[dirichletPermutation[i]];
	}
	return true;
}

void AX_HeatTransfer::updateSolution(Vector_Base<double> *x)
{
	x->store(ParametersTemperature::output->data);
}

void AX_HeatTransfer::printVersions()
{
	printParamtereStats("integration.weight", integration.weight);
	printParamtereStats("integration.N", integration.N);
	printParamtereStats("integration.dN", integration.dN);
	printParamtereStats("integration.dND", integration.dND);
	printParamtereStats("integration.jacobiDeterminant", integration.jacobiDeterminant);
	printParamtereStats("integration.jacobiInversion", integration.jacobiInversion);

	printParamtereStats("coords.node", coords.node);
	printParamtereStats("coords.gp", coords.gp);

	printParamtereStats("thickness.gp", thickness.gp);

	printParamtereStats("cooSystem.cartesian2D", cooSystem.cartesian2D);
	printParamtereStats("cooSystem.cartesian3D", cooSystem.cartesian3D);
	printParamtereStats("cooSystem.cylindric", cooSystem.cylindric);
	printParamtereStats("cooSystem.spherical", cooSystem.spherical);

	printParamtereStats("material.model.isotropic", material.model.isotropic);
	printParamtereStats("material.model.diagonal", material.model.diagonal);
	printParamtereStats("material.model.symmetric2D", material.model.symmetric2D);
	printParamtereStats("material.model.symmetric3D", material.model.symmetric3D);
	printParamtereStats("material.model.anisotropic", material.model.anisotropic);

	printParamtereStats("material.conductivityIsotropic", material.conductivityIsotropic);
	printParamtereStats("material.conductivity", material.conductivity);
	printParamtereStats("material.density", material.density);
	printParamtereStats("material.heatCapacity", material.heatCapacity);
	printParamtereStats("material.mass", material.mass);

	printParamtereStats("temp.initial.output", temp.initial.output);
	printParamtereStats("temp.initial.node", temp.initial.node);
	printParamtereStats("temp.initial.gp", temp.initial.gp);
	printParamtereStats("temp.output", temp.output);
	printParamtereStats("temp.node", temp.node);
	printParamtereStats("temp.gp", temp.gp);

	printParamtereStats("translationMotions.output", translationMotions.output);
	printParamtereStats("translationMotions.gp", translationMotions.gp);
	printParamtereStats("translationMotions.stiffness", translationMotions.stiffness);
	printParamtereStats("translationMotions.rhs", translationMotions.rhs);

	printParamtereStats("elements.stiffness", elements.stiffness);
	printParamtereStats("elements.mass", elements.mass);
	printParamtereStats("elements.rhs", elements.rhs);

	if (gradient.output)
	{
		printParamtereStats("gradient.output", gradient.output);
	}

	printParamtereStats("gradient.xi", gradient.xi);

	if (flux.output)
	{
		printParamtereStats("flux.output", flux.output);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		printf("REGION: %s\n", info::mesh->boundaryRegions[r]->name.c_str());
		printParamtereStats("convection.heatTransferCoeficient.gp", convection.heatTransferCoeficient.gp.regions[r]);
		printParamtereStats("convection.externalTemperature.gp", convection.externalTemperature.gp.regions[r]);

		printParamtereStats("heatFlow.gp", heatFlow.gp.regions[r]);
		printParamtereStats("heatFlux.gp", heatFlux.gp.regions[r]);
		printParamtereStats("q.gp", q.gp.regions[r]);
	}
}
