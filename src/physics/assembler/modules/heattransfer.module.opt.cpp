

#include "heattransfer.module.opt.h"
#include "module.opt.hpp"
#include "math/vector.dense.h"

#include "physics/assembler/operators/basis.h"
#include "physics/assembler/operators/expression.h"
#include "physics/assembler/operators/coordinates.h"
#include "physics/assembler/operators/conductivity.h"
#include "physics/assembler/operators/convection.h"
#include "physics/assembler/operators/integration.h"
#include "physics/assembler/operators/stiffness.h"
#include "physics/assembler/operators/rhs.h"
#include "physics/assembler/operators/filler.h"
#include "physics/assembler/operators/gradient.h"
#include "physics/assembler/operators/flux.h"
#include "physics/kernels/utils/geometry.h"
#include "physics/kernels/solverdataprovider/heattransfer.provider.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"

using namespace espreso;

ElementData* HeatTransferModuleOpt::phase = NULL;
ElementData* HeatTransferModuleOpt::latentHeat = NULL;

void HeatTransferModuleOpt::createParameters()
{
	if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D && info::ecf->heat_transfer_2d.kernel == HeatTransferGlobalSettings::KERNEL::OLD) {
		return;
	}
	if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D && info::ecf->heat_transfer_3d.kernel == HeatTransferGlobalSettings::KERNEL::OLD) {
		return;
	}
	ParametersTemperature::outputInitial = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
	ParametersTemperature::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
	if (info::ecf->output.results_selection.translation_motions) {
		ParametersTranslationMotions::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
	}
	if (info::ecf->output.results_selection.phase) {
		phase = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "PHASE");
	}
	if (info::ecf->output.results_selection.latent_heat) {
		latentHeat = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "LATENT_HEAT");
	}
	if (info::ecf->output.results_selection.gradient) {
		ParametersGradient::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
	}
	if (info::ecf->output.results_selection.flux) {
		ParametersFlux::output = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
	}
}

void HeatTransferModuleOpt::insertParameters(Evaluator *evaluator)
{
	for (size_t p = 0; p < evaluator->variables.size(); ++p) {
		if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", evaluator->variables[p])) {
			evaluator->params.general.push_back({ ParametersTemperature::outputInitial->data.data(), 0, 1 });
		}
		if (StringCompare::caseInsensitiveEq("TEMPERATURE", evaluator->variables[p])) {
			evaluator->params.general.push_back({ ParametersTemperature::output->data.data(), 0, 1 });
		}
		if (StringCompare::caseInsensitiveEq("COORDINATE_X", evaluator->variables[p])) {
			evaluator->params.general.push_back({ reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()), 0, 3 });
		}
		if (StringCompare::caseInsensitiveEq("COORDINATE_Y", evaluator->variables[p])) {
			evaluator->params.general.push_back({ reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()), 1, 3 });
		}
		if (info::mesh->dimension == 3) {
			if (StringCompare::caseInsensitiveEq("COORDINATE_Z", evaluator->variables[p])) {
				evaluator->params.general.push_back({ reinterpret_cast<double*>(info::mesh->nodes->coordinates->datatarray().data()), 2, 3 });
			}
		}
		if (evaluator->params.general.size() == p) {
			eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", evaluator->variables[p]);
		}
	}
}

void HeatTransferModuleOpt::initTemperature()
{
	/// This code has to solve problem that initial temperature is set to elements regions, but we need it in nodes
	/// 1. settings -> nodeInitialTemperature
	/// 2. nodeInitialTemperature -> initialTemperature (here we average values)
	/// 3. Dirichlet -> initialTemperature (if 'init_temp_respect_bc')
	/// 4. initialTemperature -> nodeInitialTemperature (TODO: check the correction with TB)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	temp.initial.node.builder = new ExpressionsToElements(temp.initial.node, 273.15);
	if (info::mesh->dimension == 2) {
		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
			temp.initial.node.setConstness(false);
		}
		examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature, *temp.initial.node.builder);
	}
	if (info::mesh->dimension == 3) {
		if (info::ecf->heat_transfer_2d.init_temp_respect_bc) {
			temp.initial.node.setConstness(false);
		}
		examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature, *temp.initial.node.builder);
	}
	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		HeatTransferModuleOpt::insertParameters(it->second.evaluator);
	}

	temp.initial.node.builder->buildAndExecute(*this);
	AverageElementsNodesToNodes(temp.initial.node, *ParametersTemperature::outputInitial).buildAndExecute(*this);
	if (info::mesh->dimension == 2 && info::ecf->heat_transfer_2d.init_temp_respect_bc) {
		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::outputInitial).buildAndExecute(*this);
	}
	if (info::mesh->dimension == 3 && info::ecf->heat_transfer_3d.init_temp_respect_bc) {
		CopyBoundaryRegionsSettingToNodes(configuration.temperature, *ParametersTemperature::outputInitial).buildAndExecute(*this);
	}
	temp.initial.node.addInputs(ParametersTemperature::output);
	temp.initial.boundary.node.addInputs(ParametersTemperature::outputInitial);
	CopyNodesToElementsNodes(*ParametersTemperature::outputInitial, temp.initial.node).now();
	CopyNodesToBoundaryNodes(*ParametersTemperature::outputInitial, temp.initial.boundary.node).now();
	ElementsGaussPointsBuilder<HeatTransferModuleOpt, 1>(integration.N, temp.initial.node, temp.initial.gp).buildAndExecute(*this);
	BoundaryGaussPointsBuilder<HeatTransferModuleOpt, 1>(integration.boundary.N, temp.initial.boundary.node, temp.initial.boundary.gp).buildAndExecute(*this);

	ParametersTemperature::output->data = ParametersTemperature::outputInitial->data;
	ParametersTemperature::output->version = ParametersTemperature::outputInitial->version;

	temp.node.addInputs(ParametersTemperature::output);
	temp.boundary.node.addInputs(ParametersTemperature::output);
	CopyElementParameters(temp.initial.node, temp.node).now();
	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node));
	builders.push_back(new ElementsGaussPointsBuilder<HeatTransferModuleOpt, 1>(integration.N, temp.node, temp.gp));
	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node));
	builders.push_back(new BoundaryGaussPointsBuilder<HeatTransferModuleOpt, 1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp));
}

HeatTransferModuleOpt::HeatTransferModuleOpt(HeatTransferModuleOpt *previous, HeatTransferLoadStepConfiguration &configuration)
: ModuleOpt(new HeatTransferSolverDataProvider(configuration)),
  configuration(configuration)
{
	geometry::computeBoundaryRegionsArea();

	double start = eslog::time();

	eslog::info("\n ============================================================================================= \n");
	eslog::info("  PHYSICS                                                                    HEAT TRANSFER 2D  \n");
	eslog::info(" ============================================================================================= \n");

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

	Basis().build(*this);
	builders.push_back(new ElementCoordinates<HeatTransferModuleOpt>(*this));
	builders.push_back(new BoundaryCoordinates<HeatTransferModuleOpt>(*this));

	builders.push_back(new ElementIntegration<HeatTransferModuleOpt>(*this));
	builders.push_back(new BoundaryIntegration<HeatTransferModuleOpt>(*this));

	if (step::loadstep == 0) {
		initTemperature(); // we need to init temperature first since other builder can read it

		if (info::mesh->dimension == 2) {
			builders.push_back(new ExpressionsToElements(thickness.gp, 1));
			examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, *thickness.gp.builder);
			builders.push_back(new ExpressionsToBoundaryFromElement(*thickness.gp.builder, thickness.boundary.gp, "THICKNESS"));

			if (configuration.translation_motions.size()) {
				switch (info::ecf->heat_transfer_2d.stabilization) {
				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
				}
				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_2d.sigma);
			}
		}
		if (info::mesh->dimension == 3) {
			thickness.gp.resize(); // dummy resize in order to avoid memory errors in some operators
			if (configuration.translation_motions.size()) {
				switch (info::ecf->heat_transfer_3d.stabilization) {
				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
				}
				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_3d.sigma);
			}
		}
		eslog::info(" ============================================================================================= \n");

		builders.push_back(new ExpressionsToElements(cooSystem.cartesian2D, 0));
		builders.push_back(new ExpressionsToElements(cooSystem.cartesian3D, 0));
		builders.push_back(new ExpressionsToElements(cooSystem.spherical, 0));
		builders.push_back(new ExpressionsToElements(cooSystem.cylindric, 0));
		builders.push_back(new ExpressionsToElements(material.model.isotropic, 1));
		builders.push_back(new ExpressionsToElements(material.model.diagonal, 1));
		builders.push_back(new ExpressionsToElements(material.model.symmetric2D, 1));
		builders.push_back(new ExpressionsToElements(material.model.symmetric3D, 1));
		builders.push_back(new ExpressionsToElements(material.model.anisotropic, 1));
		builders.push_back(new ExpressionsToElements(material.density, 1));
		builders.push_back(new ExpressionsToElements(material.heatCapacity, 1));

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
			eslog::info("                                                                                               \n");
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->heat_transfer_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->heat_transfer_3d.material_set);
		}
		builders.push_back(new ThermalConductivity(*this));
		eslog::info(" ============================================================================================= \n");
	}

	if (configuration.translation_motions.size()) {
		builders.push_back(new ExpressionsToElements(translationMotions.gp, 1));
		examineElementParameter("TRANSLATION MOTION.X", configuration.translation_motions, *translationMotions.gp.builder, 0);
		examineElementParameter("TRANSLATION MOTION.Y", configuration.translation_motions, *translationMotions.gp.builder, 1);
		if (info::mesh->dimension == 3) {
			examineElementParameter("TRANSLATION MOTION.Z", configuration.translation_motions, *translationMotions.gp.builder, 2);
		}
	}

	builders.push_back(new HeatStiffness(*this));

	if (configuration.temperature.size()) {
		builders.push_back(new ExpressionsToBoundary(dirichlet.gp));
		examineBoundaryParameter("TEMPERATURE", configuration.temperature, *dirichlet.gp.builder);
	}
	if (configuration.heat_flow.size()) {
		builders.push_back(new ExpressionsToBoundary(heatFlow.gp));
		examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, *heatFlow.gp.builder);
	}
	if (configuration.heat_flux.size()) {
		builders.push_back(new ExpressionsToBoundary(heatFlux.gp));
		examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, *heatFlux.gp.builder);
	}

	if (configuration.convection.size()) {
		builders.push_back(new ExpressionsToBoundary(convection.heatTransferCoeficient.gp));
		builders.push_back(new ExpressionsToBoundary(convection.externalTemperature.gp));

		builders.push_back(new ExpressionsToBoundary(convection.wallHeight.gp));
		builders.push_back(new ExpressionsToBoundary(convection.tiltAngle.gp));
		builders.push_back(new ExpressionsToBoundary(convection.diameter.gp));
		builders.push_back(new ExpressionsToBoundary(convection.plateLength.gp));
		builders.push_back(new ExpressionsToBoundary(convection.fluidVelocity.gp));
		builders.push_back(new ExpressionsToBoundary(convection.plateDistance.gp));
		builders.push_back(new ExpressionsToBoundary(convection.length.gp));
		builders.push_back(new ExpressionsToBoundary(convection.experimentalConstant.gp));
		builders.push_back(new ExpressionsToBoundary(convection.volumeFraction.gp));
		builders.push_back(new ExpressionsToBoundary(convection.absolutePressure.gp));

		builders.push_back(new ConvectionBuilder(*this, convection));
		examineBoundaryParameter("CONVECTION", configuration.convection, convection);
	}
//	if (configuration.diffuse_radiation.size()) {
//		eslog::info("  DIFFUSE RADIATION                                                                            \n");
//		for (auto it = configuration.diffuse_radiation.begin(); it != configuration.diffuse_radiation.end(); ++it) {
//			if (it->second.external_temperature.evaluator->parameters.size()) {
//				std::string params = Parser::join(", ", it->second.external_temperature.evaluator->parameters);
//				eslog::info("   %30s: EX. TEMP. %*s       FNC( %s )\n", it->first.c_str(), 34 - params.size(), "", params.c_str());
//			} else {
//				eslog::info("   %30s: EX. TEMP. %48g \n", it->first.c_str(), it->second.external_temperature.evaluator->eval(Evaluator::Params()));
//			}
//		}
//		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
//	}

	builders.push_back(new HeatRHS(*this));

	std::vector<OperatorBuilder*> nonempty; // there can be empty builders (e.g. from material)
	for (auto op = builders.begin(); op != builders.end(); ++op) {
		if ((*op)->build(*this)) {
			nonempty.push_back(*op);
		}
	}
	builders.swap(nonempty);

	if (gradient.output) {
		results.push_back(new Gradient(*this));
	}
	if (flux.output) {
		results.push_back(new Flux(*this));
	}

	eslog::info(" ============================================================================================= \n");
	eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	eslog::info(" ============================================================================================= \n");

	printf("ANALYZE AND BUILD: %f\n", eslog::time() - start);
}

void HeatTransferModuleOpt::nextSubstep()
{
	for (auto op = builders.begin(); op != builders.end(); ++op) {
		(*op)->now();
	}
}

void HeatTransferModuleOpt::solutionChanged(Vectors *solution)
{
	VectorDense tmp(temp.output->data.size(), temp.output->data.data());
	tmp.fillData(solution->at(0));

	temp.output->version += builders.size() + 1;
	for (auto op = builders.begin(); op != builders.end(); ++op) {
		(*op)->now();
	}
}

void HeatTransferModuleOpt::updateStiffness(double *K, esint *perm, int interval)
{
	MatricesFiller(*this, K, perm).apply(interval);
}

void HeatTransferModuleOpt::updateStiffness(double *K, esint *perm, int region, int interval)
{

}

void HeatTransferModuleOpt::updateRHS(double *RHS, esint *perm, int region, int interval)
{
	RHSFiller(*this, RHS, perm).apply(region, interval);
}

void HeatTransferModuleOpt::fillElementsInterval(int interval)
{

}

void HeatTransferModuleOpt::processSolution()
{
	for (auto op = results.begin(); op != results.end(); ++op) {
		(*op)->now();
	}
}








