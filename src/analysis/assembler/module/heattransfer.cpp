
#include "heattransfer.h"
#include "assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/conductivity.h"
#include "analysis/assembler/operators/stiffness.h"
#include "analysis/assembler/operators/filler.h"

//#include "analysis/assembler/operators/conductivity.h"
//#include "analysis/assembler/operators/convection.h"
//#include "analysis/assembler/operators/mass.h"
//#include "analysis/assembler/operators/diffusionsplit.h"
//#include "analysis/assembler/operators/advection.h"
//#include "analysis/assembler/operators/integration.h"
//#include "analysis/assembler/operators/stiffness.h"
//#include "analysis/assembler/operators/rhs.h"
//#include "analysis/assembler/operators/filler.h"
//#include "analysis/assembler/operators/gradient.h"
//#include "analysis/assembler/operators/flux.h"

#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>

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
//	if (info::ecf->output.results_selection.translation_motions) {
//		ParametersTranslationMotions::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
//	}
//	if (info::ecf->output.results_selection.phase) {
//		phase = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "PHASE");
//	}
//	if (info::ecf->output.results_selection.latent_heat) {
//		latentHeat = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "LATENT_HEAT");
//	}
//	if (info::ecf->output.results_selection.gradient) {
//		ParametersGradient::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
//	}
//	if (info::ecf->output.results_selection.flux) {
//		ParametersFlux::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
//	}

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

	Basis<AX_HeatTransfer>(*this).build();
	builders.push_back(new ElementCoordinates<AX_HeatTransfer>(*this));
	builders.push_back(new ElementIntegration<AX_HeatTransfer>(*this));

//	if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//	{
//		builders.push_back(new ElementCoordinatesSimd(*this));
//	}
//	else
//	{
//		builders.push_back(new ElementCoordinates(*this));
//	}
//
//	builders.push_back(new BoundaryCoordinates(*this));
//
//	if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//	{
//		builders.push_back(new ElementIntegrationSimd(*this));
//	}
//	else
//	{
//		builders.push_back(new ElementIntegration(*this));
//	}
//
//	builders.push_back(new BoundaryIntegration(*this));
//
//	if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//	{
//		eslog::warning("WARNING: USING EXPERIMENTAL ASSEMBLER. \n");
//		eslog::warning("WARNING: NUMBER OF ELEMENTS HAS TO BE MULTIPLE OF VECTOR LENGTH. \n");
//	}
//
	if (step::step.loadstep == 0) {
//		initTemperature(); // we need to init temperature first since other builder can read it
//
		if (info::mesh->dimension == 2) {
//			if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//			{
//				eslog::warning("WARNING: CONSTANT THICKNESS ONLY.\n");
//				builders.push_back(new ExpressionsToElementsSimd(thicknessSimd.gp, 1, "ELEMENTS THICKNESS"));
//				examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, *thicknessSimd.gp.builder);
//			}
//			// Redundant beacause of boundary
			builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, thickness.gp, 1, "ELEMENTS THICKNESS"));
			examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, *thickness.gp.builder);
//			builders.push_back(new ExpressionsToBoundaryFromElement(*thickness.gp.builder, thickness.boundary.gp, "BOUDARY THICKNESS"));
//
//			if (configuration.translation_motions.size()) {
//				switch (info::ecf->heat_transfer_2d.stabilization) {
//				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
//				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
//				}
//				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_2d.sigma);
//			}
		}
//		if (info::mesh->dimension == 3) {
//			thickness.gp.resize(); // dummy resize in order to avoid memory errors in some operators
//			if (configuration.translation_motions.size()) {
//				switch (info::ecf->heat_transfer_3d.stabilization) {
//				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
//				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
//				}
//				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_3d.sigma);
//			}
//		}
//		eslog::info(" ============================================================================================= \n");
//
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			builders.push_back(new ExpressionsToElementsSimd(materialSimd.model.isotropic, 1, "ISOTROPIC MATERIAL MODEL"));
//		}

		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, cooSystem.cartesian2D, 0, "CARTESIAN 2D COORDINATE SYSTEM"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, cooSystem.cartesian3D, 0, "CARTESIAN 3D COORDINATE SYSTEM"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, cooSystem.spherical, 0, "SPHERICAL COORDINATE SYSTEM"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, cooSystem.cylindric, 0, "CYLINDRIC COORDINATE SYSTEM"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.model.isotropic, 1, "ISOTROPIC MATERIAL MODEL"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.model.diagonal, 1, "DIAGONAL MATERIAL MODEL"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.model.symmetric2D, 1, "SYMMETRIC MATERIAL MODEL"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.model.symmetric3D, 1, "SYMMETRIC MATERIAL MODEL"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.model.anisotropic, 1, "ANISOTROPIC MATERIAL MODEL"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.density, 1, "MATERIAL DENSITY"));
		builders.push_back(new ExpressionsToElementsSimple<AX_HeatTransfer>(*this, material.heatCapacity, 1, "MATERIAL HEAT CAPACITY"));

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
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, *cooSystem.cartesian2D.builder, 0);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x, *cooSystem.cartesian3D.builder, 0);
					examineMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y, *cooSystem.cartesian3D.builder, 1);
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, *cooSystem.cartesian3D.builder, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, *cooSystem.spherical.builder, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, *cooSystem.spherical.builder, 1);
					examineMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z, *cooSystem.spherical.builder, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, *cooSystem.cylindric.builder, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, *cooSystem.cylindric.builder, 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, *cooSystem.cylindric.builder, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, *cooSystem.cylindric.builder, 1);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			examineMaterialParameter(mat->name, "DENSITY", mat->density, *material.density.builder, 0);
			examineMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity, *material.heatCapacity.builder, 0);
			eslog::info("                                                                                               \n");

			if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
			{
				eslog::warning("WARNING: CONSTANT THERMAL CONDUCTIVITY ONLY.\n");
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
					eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *materialSimd.model.isotropic.builder, 0);
					break;
				default:
					eslog::error("EXPERIMENTAL ASSEMBLER SUPPORTS ISOTROPIC MODEL OF THERMAL CONDCTIVITY ONLY.\n");
				}
			}
			else
			{
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
					eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.isotropic.builder, 0);
					break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:
					eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
					if (info::mesh->dimension == 2) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.diagonal.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.diagonal.builder, 1);
					}
					if (info::mesh->dimension == 3) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.diagonal.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.diagonal.builder, 1);
						examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), *material.model.diagonal.builder, 2);
					}
					break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
					eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
					if (info::mesh->dimension == 2) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.symmetric2D.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.symmetric2D.builder, 1);
						examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), *material.model.symmetric2D.builder, 2);
					}
					if (info::mesh->dimension == 3) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.symmetric3D.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.symmetric3D.builder, 1);
						examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), *material.model.symmetric3D.builder, 2);
						examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), *material.model.symmetric3D.builder, 3);
						examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), *material.model.symmetric3D.builder, 4);
						examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), *material.model.symmetric3D.builder, 5);
					}
					break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
					eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
					if (info::mesh->dimension == 2) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.anisotropic.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.anisotropic.builder, 1);
						examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), *material.model.anisotropic.builder, 2);
						examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0), *material.model.anisotropic.builder, 3);
					}
					if (info::mesh->dimension == 3) {
						examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), *material.model.anisotropic.builder, 0);
						examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), *material.model.anisotropic.builder, 1);
						examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), *material.model.anisotropic.builder, 2);
						examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), *material.model.anisotropic.builder, 3);
						examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), *material.model.anisotropic.builder, 4);
						examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), *material.model.anisotropic.builder, 5);
						examineMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0), *material.model.anisotropic.builder, 6);
						examineMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1), *material.model.anisotropic.builder, 7);
						examineMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0), *material.model.anisotropic.builder, 8);
					}
					break;
				}
			}
			eslog::info("                                                                                               \n");
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->heat_transfer_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->heat_transfer_3d.material_set);
		}
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			builders.push_back(new ThermalConductivitySimd(*this));
//		}
//		else
//		{
//			builders.push_back(new ThermalConductivity(*this));
//		}
		builders.push_back(new ThermalConductivity<AX_HeatTransfer>(*this));

		eslog::info(" ============================================================================================= \n");
	}

	gradient.xi.resize(1);

//
//	builders.push_back(new MaterialMassBuilder(*this));
//
//	if (info::mesh->dimension == 2 && info::ecf->heat_transfer_2d.diffusion_split) {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			eslog::error("EXPERIMENTAL ASSEMBLER DOESNT SUPPORT GRADIENT.\n");
//		}
//		Gradient(*this).buildAndExecute(*this);
//		builders.push_back(new DiffusionSplit(*this));
//	} else {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			gradientSimd.xi.resizeAligned(SIMD::size*sizeof(double), 1);
//		}
//		else
//		{
//			gradient.xi.resize(1);
//		}
//	}
//	if (info::mesh->dimension == 3 && info::ecf->heat_transfer_3d.diffusion_split) {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			eslog::error("EXPERIMENTAL ASSEMBLER DOESNT SUPPORT GRADIENT.\n");
//		}
//		Gradient(*this).buildAndExecute(*this);
//		builders.push_back(new DiffusionSplit(*this));
//	} else {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			gradientSimd.xi.resizeAligned(SIMD::size*sizeof(double), 1);
//		}
//		else
//		{
//			gradient.xi.resize(1);
//		}
//	}
//
//	if (configuration.translation_motions.size()) {
//		builders.push_back(new ExpressionsToElementsSimple(translationMotions.gp, 1, "TRANSLATION MOTION"));
//		examineElementParameter("TRANSLATION MOTION.X", configuration.translation_motions, *translationMotions.gp.builder, 0);
//		examineElementParameter("TRANSLATION MOTION.Y", configuration.translation_motions, *translationMotions.gp.builder, 1);
//		if (info::mesh->dimension == 3) {
//			examineElementParameter("TRANSLATION MOTION.Z", configuration.translation_motions, *translationMotions.gp.builder, 2);
//		}
//		builders.push_back(new TranslationMotion(*this));
//	}
//
//	if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//	{
//		builders.push_back(new HeatStiffnessSimd(*this));
//	}
//	else
//	{
//		builders.push_back(new HeatStiffness(*this));
//	}

	builders.push_back(new HeatStiffness<AX_HeatTransfer>(*this));

//
//	if (configuration.temperature.size()) {
//		builders.push_back(new ExpressionsToBoundary(dirichlet.gp, "BOUNDARY TEMPERATURE"));
//		examineBoundaryParameter("TEMPERATURE", configuration.temperature, *dirichlet.gp.builder);
//	}
//	if (configuration.heat_flow.size()) {
//		builders.push_back(new ExpressionsToBoundary(heatFlow.gp, "BOUNDARY HEAT FLOW"));
//		examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, *heatFlow.gp.builder);
//	}
//	if (configuration.heat_flux.size()) {
//		builders.push_back(new ExpressionsToBoundary(heatFlux.gp, "BOUNDARY HEAT FLUX"));
//		examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, *heatFlux.gp.builder);
//	}
//
//	if (configuration.convection.size()) {
//		builders.push_back(new ExpressionsToBoundary(convection.heatTransferCoeficient.gp, "CONVECTION: HTC"));
//		builders.push_back(new ExpressionsToBoundary(convection.externalTemperature.gp, "CONVECTION: EXTERNAL TEMPERATURE"));
//
//		builders.push_back(new ExpressionsToBoundary(convection.wallHeight.gp, "CONVECTION: WALL HEIGHT"));
//		builders.push_back(new ExpressionsToBoundary(convection.tiltAngle.gp, "CONVECTION: TILT ANGLE"));
//		builders.push_back(new ExpressionsToBoundary(convection.diameter.gp, "CONVECTION: DIAMETER"));
//		builders.push_back(new ExpressionsToBoundary(convection.plateLength.gp, "CONVECTION: PLATE LENGTH"));
//		builders.push_back(new ExpressionsToBoundary(convection.fluidVelocity.gp, "CONVECTION: FLUID VELOCITY"));
//		builders.push_back(new ExpressionsToBoundary(convection.plateDistance.gp, "CONVECTION: PLATE DISTANCE"));
//		builders.push_back(new ExpressionsToBoundary(convection.length.gp, "CONVECTION: LENGTH"));
//		builders.push_back(new ExpressionsToBoundary(convection.experimentalConstant.gp, "CONVECTION: EXPERIMENTAL CONSTANT"));
//		builders.push_back(new ExpressionsToBoundary(convection.volumeFraction.gp, "CONVECTION: VOLUME FRACTION"));
//		builders.push_back(new ExpressionsToBoundary(convection.absolutePressure.gp, "CONVECTION: ABSOLUTE PRESSURE"));
//
//		builders.push_back(new ConvectionBuilder(*this, convection));
//		examineBoundaryParameter("CONVECTION", configuration.convection, convection);
//	}
////	if (configuration.diffuse_radiation.size()) {
////		eslog::info("  DIFFUSE RADIATION                                                                            \n");
////		for (auto it = configuration.diffuse_radiation.begin(); it != configuration.diffuse_radiation.end(); ++it) {
////			if (it->second.external_temperature.evaluator->parameters.size()) {
////				std::string params = Parser::join(", ", it->second.external_temperature.evaluator->parameters);
////				eslog::info("   %30s: EX. TEMP. %*s       FNC( %s )\n", it->first.c_str(), 34 - params.size(), "", params.c_str());
////			} else {
////				eslog::info("   %30s: EX. TEMP. %48g \n", it->first.c_str(), it->second.external_temperature.evaluator->eval(Evaluator::Params()));
////			}
////		}
////		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
////	}
//
//	builders.push_back(new HeatRHS(*this));
//
	std::vector<OperatorBuilder*> nonempty; // there can be empty builders (e.g. from material)
	for (auto op = builders.begin(); op != builders.end(); ++op) {
		if ((*op)->build()) {
			nonempty.push_back(*op);
		} else {
			if (Operator::print > 1) printf("DROP %s\n", (*op)->name);
		}
	}
	builders.swap(nonempty);
//
//	if (gradient.output) {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			eslog::error("EXPERIMENTAL ASSEMBLER DOESNT SUPPORT GRADIENT OUTPUT.\n");
//		}
//		results.push_back(new Gradient(*this));
//	}
//	if (flux.output) {
//		if (gsettings.kernel == HeatTransferGlobalSettings::KERNEL::VEC)
//		{
//			eslog::error("EXPERIMENTAL ASSEMBLER DOESNT SUPPORT FLUX OUTPUT.\n");
//		}
//		results.push_back(new Flux(*this));
//	}

	eslog::info(" ============================================================================================= \n");
	eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	eslog::info(" ============================================================================================= \n");
}

void AX_HeatTransfer::next()
{
	updateVersions();
	for (auto op = builders.begin(); op != builders.end(); ++op) {
		(*op)->now();
	}
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
		MatricesFiller<AX_HeatTransfer>(*this, mapping.elements[i].data, mapping.elements[i].position).apply(i);
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
