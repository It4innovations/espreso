
#include "heattransfer.h"
#include "heattransfer.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/temperature.h"
#include "analysis/assembler/operators/advection.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/heattransfer.f.h"
#include "analysis/assembler/operators/heattransfer.K.h"
#include "analysis/assembler/operators/filler.h"
#include "analysis/assembler/operators/gradient.h"
#include "analysis/assembler/operators/flux.h"

#include "basis/expression/variable.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"

#include <numeric>
#include <algorithm>

namespace espreso {

NodeData* HeatTransfer::Results::temperature = nullptr;
NodeData* HeatTransfer::Results::initialTemperature = nullptr;
ElementData* HeatTransfer::Results::translationMotion = nullptr;
ElementData* HeatTransfer::Results::gradient = nullptr;
ElementData* HeatTransfer::Results::flux = nullptr;

HeatTransfer::HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
	elements.stiffness.setConstness(false);
	elements.stiffness.resize();
	elements.rhs.setConstness(false);
	elements.rhs.resize();
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		elements.boundary.rhs.regions[r].setConstness(false);
		elements.boundary.rhs.regions[r].resize();
		if (configuration.heat_flow.end() != configuration.heat_flow.find(info::mesh->boundaryRegions[r]->name)) {
			bfilter[r] = 1;
		}
		if (configuration.heat_flux.end() != configuration.heat_flux.find(info::mesh->boundaryRegions[r]->name)) {
			bfilter[r] = 1;
		}
	}
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
		correct &= checkBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature);
		generateBoundaryExpression<ExternalNodeExpression>(boundaryOps, configuration.temperature, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.temperature[n][s] = value; });
	}

//	if (settings.init_temp_respect_bc) {
//		temp.initial.node.setConstness(false);
//	}
//	examineElementParameter("INITIAL TEMPERATURE", settings.initial_temperature, temp.initial.node.externalValues);
//	// TODO: fix
////	fromExpression(*this, temp.initial.node, temp.initial.node.externalValues);
//	_evaluate();
//
////	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
////		HeatTransfer::insertParameters(it->second.evaluator);
////	}
////
////	temp.initial.node.builder->buildAndExecute(*this);
//
//	// TODO: fix
////	averageEnodesToNodes(temp.initial.node, *Results::initialTemperature);
////	Results::temperature->data = Results::initialTemperature->data;
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
//	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end() || Results::gradient) {
//		copyNodesToEnodes(*this, *Results::temperature, temp.node);
//	}
//
//	if (Variable::list.egps.find("TEMPERATURE") != Variable::list.egps.end()) {
//		moveEnodesToGPs(*this, temp.node, temp.gp, 1);
//		Variable::list.egps["TEMPERATURE"] = new ParameterVariable(temp.gp.data, temp.gp.isconst, temp.gp.update, 0, 1);
//	}
//
//
//
////	ParametersTemperature::output->data = ParametersTemperature::Initial::output->data;
////	CopyElementParameters(temp.initial.node, temp.node, "COPY INITIAL TEMPERATURE TO ELEMENT NODES").buildAndExecute(*this);
////	builders.push_back(new CopyNodesToElementsNodes(*ParametersTemperature::output, temp.node, "COPY TEMPERATURE TO ELEMENTS NODES"));
////	builders.push_back(new ElementsGaussPointsBuilder<1>(integration.N, temp.node, temp.gp, "INTEGRATE TEMPERATURE INTO ELEMENTS GAUSS POINTS"));
////	builders.push_back(new CopyNodesToBoundaryNodes(*ParametersTemperature::output, temp.boundary.node, "COPY TEMPERATURE TO BOUNDARY NODES"));
////	builders.push_back(new BoundaryGaussPointsBuilder<1>(integration.boundary.N, temp.boundary.node, temp.boundary.gp, "INTEGRATE TEMPERATURE INTO BOUNDARY GAUSS POINTS"));
//
//	results();
	return correct;
}

void HeatTransfer::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);

	initParameters();

	eslog::info(" ============================================================================================= \n");
	bool correct = true;
	correct &= initTemperature();

	if (step::step.loadstep == 0) {
		if (info::mesh->dimension == 2) {
			correct &= checkElementParameter("THICKNESS", settings.thickness);
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
					correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("ROTATION.X", mat->coordinate_system.rotation.x);
					correct &= checkExpression("ROTATION.Y", mat->coordinate_system.rotation.y);
					correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
					correct &= checkExpression("CENTER.Z", mat->coordinate_system.center.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= checkExpression("DENSITY", mat->density);
			correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(1, 0));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
					correct &= checkExpression("KYX", mat->thermal_conductivity.values.get(1, 0));
					correct &= checkExpression("KZY", mat->thermal_conductivity.values.get(2, 1));
					correct &= checkExpression("KZX", mat->thermal_conductivity.values.get(2, 0));
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		if (configuration.translation_motions.size()) {
			for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
				const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:   etype[interval] = HeatTransferElementType::ASYMMETRIC_ISOTROPIC; break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:    etype[interval] = HeatTransferElementType::ASYMMETRIC_GENERAL  ; break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   etype[interval] = HeatTransferElementType::ASYMMETRIC_GENERAL  ; break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: etype[interval] = HeatTransferElementType::ASYMMETRIC_GENERAL  ; break;
				}
			}
		} else {
			for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
				const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:   etype[interval] = HeatTransferElementType::SYMMETRIC_ISOTROPIC ; break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:    etype[interval] = HeatTransferElementType::SYMMETRIC_GENERAL   ; break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   etype[interval] = HeatTransferElementType::SYMMETRIC_GENERAL   ; break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: etype[interval] = HeatTransferElementType::ASYMMETRIC_GENERAL  ; break;
				}
			}
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			switch (info::mesh->boundaryRegions[r]->dimension) {
			case 0: std::fill(btype[r].begin(), btype[r].end(), HeatTransferElementType::NODE); break;
			case 1: std::fill(btype[r].begin(), btype[r].end(), HeatTransferElementType::EDGE); break;
			case 2: std::fill(btype[r].begin(), btype[r].end(), HeatTransferElementType::FACE); break;
			}
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);
		eslog::info(" ============================================================================================= \n");
	}

	// we cannot generate functions since 'etype' is filled
	generateBaseFunctions(etype, elementOps);
	generateBaseFunctions(bfilter, boundaryOps);

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin;
		bool cooToGPs = false;

		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			cooToGPs |= mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
		}

		if (cooToGPs) {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodesAndGPs>(i, etype[i], procNodes));
		} else {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodes>(i, etype[i], procNodes));
		}

		if (Results::gradient != nullptr || Results::flux != nullptr) {
			elementOps[i].push_back(generateElementOperator<TemperatureToElementNodes>(i, etype[i], procNodes, Results::temperature->data.data()));
		}
	}

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
					boundaryOps[r][i].push_back(generateBoundaryOperator<CoordinatesToElementNodes>(r, i, procNodes));
				}
			}
		}
	}

	if (info::mesh->dimension == 2) {
		generateElementExpression2D<ExternalGPsExpression>(etype, elementOps, settings.thickness, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.thickness[gp][s] = value; });
	}
	generateConductivity();
	generateElementOperators<Integration>(etype, elementOps);
	generateBoundaryOperators<Integration>(bfilter, boundaryOps);
	volume();

	if (configuration.translation_motions.size()) { // it updates conductivity
		correct &= checkElementParameter("TRANSLATION MOTIONS", configuration.translation_motions);
		generateElementAsymmetricTypeExpression<ExternalGPsExpression>(etype, elementOps, configuration.translation_motions, 0, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][0][s] = value; });
		generateElementAsymmetricTypeExpression<ExternalGPsExpression>(etype, elementOps, configuration.translation_motions, 1, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][1][s] = value; });
		if (info::mesh->dimension == 3) {
			generateElementAsymmetricTypeExpression<ExternalGPsExpression>(etype, elementOps, configuration.translation_motions, 2, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][2][s] = value; });
		}
		generateElementAsymmetricOperators<Advection>(etype, elementOps, elements.stiffness);
	}

//	gradient.xi.resize(1);
//	controller.prepare(gradient.xi);
	generateElementOperators<HeatTransferStiffness>(etype, elementOps, elements.stiffness);
	if (configuration.heat_source.size()) {
		correct &= checkElementParameter("HEAT SOURCE", configuration.heat_source);
		generateElementExpression<ExternalGPsExpression>(etype, elementOps, configuration.heat_source, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatSource[gp][s] = value; });
		generateElementOperators<HeatSource>(etype, elementOps, elements.rhs);
	}
	if (configuration.heat_flow.size()) {
		correct &= checkBoundaryParameter("HEAT FLOW", configuration.heat_flow);
		generateBoundaryExpression<ExternalGPsExpression>(boundaryOps, configuration.heat_flow, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatFlow[gp][s] = value; });
	}
	if (configuration.heat_flux.size()) {
		correct &= checkBoundaryParameter("HEAT FLUX", configuration.heat_flux);
		generateBoundaryExpression<ExternalGPsExpression>(boundaryOps, configuration.heat_flux, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatFlux[gp][s] = value; });
	}
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				boundaryOps[r][i].push_back(generateBoundaryOperator<BoundaryHeat>(r, i, elements.boundary.rhs.regions[r]));
			}
		}
	}

	if (Results::gradient != nullptr) {
		generateElementOperators<TemperatureGradient>(etype, elementOps, Results::gradient);
	}
	if (Results::flux != nullptr) {
		generateElementOperators<TemperatureFlux>(etype, elementOps, Results::flux);
	}

	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
	eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize());
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
	this->K = scheme.K;
	switch (scheme.K->shape) {
	case Matrix_Shape::FULL:  generateElementOperators<GeneralMatricFiller>(etype, elementOps, 1, elements.stiffness, scheme.K); break;
	case Matrix_Shape::UPPER: generateElementOperators<SymmetricMatricFiller>(etype, elementOps, 1, elements.stiffness, scheme.K); break;
	}

	generateElementOperators<VectorFiller>(etype, elementOps, 1, elements.rhs, scheme.f);

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (configuration.temperature.end() == configuration.temperature.find(info::mesh->boundaryRegions[r]->name)) {
			if (bfilter[r]) {
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 0:
					for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
						boundaryOps[r][t].push_back(generateNodeFiller<VectorFiller>(r, t, 1, elements.rhs, scheme.f));
					}
					break;
				case 1:
					for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
						boundaryOps[r][i].push_back(generateEdgeFiller<VectorFiller>(r, i, 1, elements.boundary.rhs.regions[r], scheme.f));
					}
					break;
				case 2:
					for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
						boundaryOps[r][i].push_back(generateFaceFiller<VectorFiller>(r, i, 1, elements.boundary.rhs.regions[r], scheme.f));
					}
					break;
				}
			}
		} else {
			// DIRICHLET
			for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				boundaryOps[r][t].push_back(generateNodeSetter<VectorSetter>(r, t, 1, scheme.dirichlet, [] (auto &element, const size_t &n, const size_t &dof, const size_t &s) { return element.temperature[n][s]; }));
			}
		}
	}
}

void HeatTransfer::evaluate(SteadyState &scheme)
{
	reset(scheme.K, scheme.f, scheme.dirichlet);
	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n", assemble(ActionOperator::Action::ASSEMBLE));
	std::fill(elements.stiffness.data->datatarray().begin(), elements.stiffness.data->datatarray().end(), 0);
	eslog::info("       = SIMD LOOP REASSEMBLE                                           %12.8f s = \n", assemble(ActionOperator::Action::REASSEMBLE));
	eslog::info("       = FILL MATRICES                                                  %12.8f s = \n", assemble(ActionOperator::Action::FILL));
	update(scheme.K, scheme.f);
}

void HeatTransfer::volume()
{
	std::vector<double> evolume(info::mesh->elements->eintervals.size());
	std::vector<double> bvolume(info::mesh->boundaryRegions.size());

	generateElementOperators<Volume>(etype, elementOps, evolume);
	generateBoundaryOperators<Volume>(bfilter, boundaryOps, bvolume);
	assemble(ActionOperator::Action::ASSEMBLE);
	dropLastOperators(elementOps);
	dropLastOperators(bfilter, boundaryOps);

	printElementVolume(evolume);
	printBoundarySurface(bvolume);
}

size_t HeatTransfer::esize()
{
	size_t max = 0;
	std::vector<size_t> esize(info::mesh->elements->eintervals.size());

	generateElementOperators<DataDescriptorElementSize>(etype, elementOps, esize);
	for (size_t i = 0; i < elementOps.size(); ++i) {
		max = std::max(max, esize[i]);
	}
	dropLastOperators(elementOps);
	return max;
}

void HeatTransfer::updateSolution(SteadyState &scheme)
{
	scheme.x->storeTo(Results::temperature->data);
	assemble(ActionOperator::Action::SOLUTION);
//	temp.node.setUpdate(1);
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
double HeatTransfer::operatorsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	if (this->K == nullptr) {
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
	if (elements == 0) return 0;
	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->simd(element);
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
	CoordinatesToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > coo(interval, procNodes);
	TemperatureToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > temp(interval, procNodes, Results::temperature->data.data());
	Integration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > integration(interval);
	HeatTransferStiffness<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > stiffness(interval, this->elements.stiffness);

	SymmetricMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > symFiller(interval, 1, this->elements.stiffness, this->K);
	GeneralMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > asymFiller(interval, 1, this->elements.stiffness, this->K);

	TemperatureGradient<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > gradient(interval, Results::gradient);
	TemperatureFlux<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > flux(interval, Results::flux);

	coo.move(SIMD::size);
	temp.move(SIMD::size);
	stiffness.move(SIMD::size);
	symFiller.move(SIMD::size);

	gradient.move(ndim * SIMD::size);
	flux.move(ndim * SIMD::size);

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool constConductivity = true, rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = info::ecf->output.results_selection.gradient;
	bool computeFlux = info::ecf->output.results_selection.flux;
	bool getTemp = action == ActionOperator::SOLUTION && (computeGradient || computeFlux);
	bool isSymmetric = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ANISOTROPIC && (configuration.translation_motions.size() == 0 || settings.sigma == 0);

	double start = eslog::time();
	esint chunks = elements / SIMD::size;
	for (esint c = 1; c < chunks; ++c) {
		coo.simd(element);
		integration.simd(element);
		if (getTemp) {
			temp.simd(element);
		}
		if (!constConductivity) {

		}
		if (!rotateConductivity) {

		}
		if (computeK) {
			stiffness.simd(element);
		}
		if (action == ActionOperator::FILL) {
			if (isSymmetric) {
				symFiller.simd(element);
			} else {
				asymFiller.simd(element);
			}
		}
		if (computeGradient) {
			gradient.simd(element);
		}
		if (computeFlux) {
			flux.simd(element);
		}
	}
	double end = eslog::time();

	if (elements % SIMD::size) {
		eslog::error("peel loop is not supported\n");
		// peel is never needed
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-(int)std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-(esint)SIMD::size);
			}
		}
	}
	return end - start;
}

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
double HeatTransfer::manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	eslog::info("       = LOOP TYPE                                                           MANUAL = \n");
	if (elements == 0) return 0;
	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;
	std::vector<DataDescriptor<nodes, gps, ndim, edim, etype>*> active; active.reserve(ops.size());

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				if ((*op)->isconst) {
					dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->simd(element);
				} else {
					active.push_back(dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op));
					active.back()->simd(element);
				}
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->peel(element, elements);
			}
		}
	}

	double start = eslog::time();
	esint chunks = elements / SIMD::size;
	for (esint c = 1; c < chunks; ++c) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->simd(element);
		}
	}
	double end = eslog::time();

	if (elements % SIMD::size) {
		for (auto op = active.cbegin(); op != active.cend(); ++op) {
			(*op)->peel(element, elements % SIMD::size);
		}
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-(int)std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-elements);
			}
		}
	}
	return end - start;
}

template <>
double HeatTransfer::instantiate2D<HeatTransferElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<HeatTransferDataDescriptor, 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>(action, ops, elements);
}

template <>
double HeatTransfer::instantiate2D<HeatTransferElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<HeatTransferDataDescriptor, 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	default: return 0;
	}
}

template <int etype>
double HeatTransfer::instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (settings.loop) {
	case PhysicsConfiguration::LOOP::INHERITANCE:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, elements); break;
		default: return 0;
		};
	case PhysicsConfiguration::LOOP::OPERATORS:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TRIANGLE3): return operatorsloop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::TRIANGLE6): return operatorsloop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE4):   return operatorsloop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE8):   return operatorsloop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
		default: return 0;
		};
	case PhysicsConfiguration::LOOP::MANUAL:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TRIANGLE3): return manualloop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::TRIANGLE6): return manualloop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE4):   return manualloop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE8):   return manualloop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
		default: return 0;
		};
	}
	return 0;
}

template <>
double HeatTransfer::instantiate3D<HeatTransferElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<HeatTransferDataDescriptor, 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>(action, ops, elements);
}

template <>
double HeatTransfer::instantiate3D<HeatTransferElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<HeatTransferDataDescriptor, 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	default: return 0;
	}
}

template <>
double HeatTransfer::instantiate3D<HeatTransferElementType::FACE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	default: return 0;
	}
}

template <int etype>
double HeatTransfer::instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (settings.loop) {
	case PhysicsConfiguration::LOOP::INHERITANCE:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TETRA4):    return loop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::TETRA10):   return loop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID5):  return loop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID13): return loop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA6):   return loop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA15):  return loop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA8):     return loop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA20):    return loop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, elements); break;
		default: return 0;
		};
	case PhysicsConfiguration::LOOP::OPERATORS:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TETRA4):    return operatorsloop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::TETRA10):   return operatorsloop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID5):  return operatorsloop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID13): return operatorsloop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA6):   return operatorsloop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA15):  return operatorsloop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA8):     return operatorsloop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA20):    return operatorsloop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
		default: return 0;
		};
	case PhysicsConfiguration::LOOP::MANUAL:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TETRA4):    return manualloop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::TETRA10):   return manualloop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID5):  return manualloop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID13): return manualloop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA6):   return manualloop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA15):  return manualloop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA8):     return manualloop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA20):    return manualloop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
		default: return 0;
		};
	}
	return 0;
}

double HeatTransfer::instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiate2D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiate2D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiate2D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiate2D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);

		// boundary
		case HeatTransferElementType::EDGE: return instantiate2D<HeatTransferElementType::EDGE>(action, code, ops, interval, elements);
		case HeatTransferElementType::NODE: return instantiate2D<HeatTransferElementType::NODE>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiate3D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiate3D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiate3D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiate3D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);

		// boundary
		case HeatTransferElementType::FACE: return instantiate3D<HeatTransferElementType::FACE>(action, code, ops, interval, elements);
		case HeatTransferElementType::EDGE: return instantiate3D<HeatTransferElementType::EDGE>(action, code, ops, interval, elements);
		case HeatTransferElementType::NODE: return instantiate3D<HeatTransferElementType::NODE>(action, code, ops, interval, elements);
		}
	}
	return 0;
}

}

