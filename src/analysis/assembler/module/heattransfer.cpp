
#include "heattransfer.h"
#include "heattransfer.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/temperature.h"
#include "analysis/assembler/operators/initialtemperature.h"
#include "analysis/assembler/operators/thickness.h"
#include "analysis/assembler/operators/advection.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/conductivity.coordinatesystem.h"
#include "analysis/assembler/operators/heattransfer.f.h"
#include "analysis/assembler/operators/heattransfer.K.h"
#include "analysis/assembler/operators/filler.h"
#include "analysis/assembler/operators/gradient.h"
#include "analysis/assembler/operators/flux.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"
#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

#include <iostream>
namespace espreso {

NodeData* HeatTransfer::Results::temperature = nullptr;
NodeData* HeatTransfer::Results::initialTemperature = nullptr;
NodeData* HeatTransfer::Results::thickness = nullptr;
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
		if (configuration.temperature.end() != configuration.temperature.find(info::mesh->boundaryRegions[r]->name)) {
			bfilter[r] = 1;
		}
	}
	cossin_conditions.resize(info::mesh->elements->eintervals.size());
}

void HeatTransfer::initParameters()
{
	if (Results::temperature == nullptr) {
		Results::temperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
	}
	if (Results::initialTemperature == nullptr) {
		Results::initialTemperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
	}
	if (Results::thickness == nullptr && info::mesh->dimension == 2) {
		Results::thickness = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "THICKNESS");
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

	if (configuration.temperature.size()) {
		correct &= checkBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature);
	}

	if (settings.initial_temperature.size()) {
		correct &= checkElementParameter("INITIAL TEMPERATURE", settings.initial_temperature);
	}

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
	}

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
					boundaryOps[r][i].push_back(generateBoundaryOperator<CoordinatesToElementNodes>(r, i, procNodes));
				}
			} else {
				for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					auto procNodes = info::mesh->boundaryRegions[r]->nodes->cbegin(t);
					boundaryOps[r][t].push_back(generateBoundaryNodeOperator<CoordinatesToElementNodes>(r, t, procNodes));
				}
			}
		}
	}

	if (info::mesh->dimension == 2) {
		generateElementExpression2D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, settings.thickness, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.thickness[gp][s] = value; });
	}

	initTemperatureAndThickness();

	if (info::mesh->dimension == 2) {
		for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (bfilter[r]) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
						auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
						boundaryOps[r][i].push_back(generateBoundaryEdge2DOperator<ThicknessToElementNodes>(r, i, procNodes, Results::thickness->data.data()));
					}
				} else {
					for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
						auto procNodes = info::mesh->boundaryRegions[r]->nodes->cbegin(t);
						boundaryOps[r][t].push_back(generateBoundaryNode2DOperator<ThicknessToElementNodes>(r, t, procNodes, Results::thickness->data.data()));
					}
				}
			}
		}
	}

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin;
		bool tempToGPs = false;

		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		for (size_t v = 0; !tempToGPs && v < mat->thermal_conductivity.values.size * mat->thermal_conductivity.values.size; ++v) {
			if (mat->thermal_conductivity.values.values[v].evaluator) {
				tempToGPs |= mat->thermal_conductivity.values.values[v].evaluator->needTemperature(info::mesh->elements->eintervals[i].thread);
			}
		}
		if (tempToGPs) {
			elementOps[i].push_back(generateElementOperator<TemperatureToElementNodesAndGPs>(i, etype[i], procNodes, Results::temperature->data.data()));
		} else {
			if (Results::gradient != nullptr || Results::flux != nullptr) {
				elementOps[i].push_back(generateElementOperator<TemperatureToElementNodes>(i, etype[i], procNodes, Results::temperature->data.data()));
				elementOps[i].back()->action = ActionOperator::SOLUTION;
			}
		}
	}
	generateConductivity();
	generateElementOperators<Integration>(etype, elementOps);
	generateBoundaryOperators<Integration>(bfilter, boundaryOps);
	volume();

	if (configuration.temperature.size()) {
		generateBoundaryExpression<ExternalNodeExpressionWithCoordinates>(boundaryOps, configuration.temperature, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.temp[n][s] = value; });
	}

	if (configuration.translation_motions.size()) { // it updates conductivity
		correct &= checkElementParameter("TRANSLATION MOTIONS", configuration.translation_motions);
		generateElementAsymmetricTypeExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.translation_motions, 0, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][0][s] = value; });
		generateElementAsymmetricTypeExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.translation_motions, 1, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][1][s] = value; });
		if (info::mesh->dimension == 3) {
			generateElementAsymmetricTypeExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.translation_motions, 2, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.advection[gp][2][s] = value; });
		}
		generateElementAsymmetricOperators<Advection>(etype, elementOps, elements.stiffness);
	}

//	gradient.xi.resize(1);
//	controller.prepare(gradient.xi);
	generateElementOperators<HeatTransferStiffness>(etype, elementOps, elements.stiffness);
	if (configuration.heat_source.size()) {
		correct &= checkElementParameter("HEAT SOURCE", configuration.heat_source);
		generateElementExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.heat_source, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatSource[gp][s] = value; });
		generateElementOperators<HeatSource>(etype, elementOps, elements.rhs);
	}
	if (configuration.heat_flow.size()) {
		correct &= checkBoundaryParameter("HEAT FLOW", configuration.heat_flow);
		generateBoundaryExpression<ExternalGpsExpressionWithCoordinates>(boundaryOps, configuration.heat_flow, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatFlow[gp][s] = value; });
	}
	if (configuration.heat_flux.size()) {
		correct &= checkBoundaryParameter("HEAT FLUX", configuration.heat_flux);
		generateBoundaryExpression<ExternalGpsExpressionWithCoordinates>(boundaryOps, configuration.heat_flux, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatFlux[gp][s] = value; });
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

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		if (Results::gradient == nullptr && Results::flux == nullptr) {
			for (size_t j = 0; j < elementOps[i].size(); ++j) {
				ActionOperator::removeSolution(elementOps[i][j]->action);
			}
		}
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
	this->f = scheme.f;
	switch (scheme.K->shape) {
	case Matrix_Shape::FULL:
		for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			switch (etype[i]) {
			case HeatTransferElementType::SYMMETRIC_ISOTROPIC:
			case HeatTransferElementType::SYMMETRIC_GENERAL:
				elementOps[i].push_back(generateElementOperator<SymmetricToFullMatricFiller>(i, etype[i], 1, elements.stiffness, scheme.K)); break;
				break;
			default:
				elementOps[i].push_back(generateElementOperator<GeneralMatricFiller>(i, etype[i], 1, elements.stiffness, scheme.K)); break;
				break;
			}
		}
		break;
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
				boundaryOps[r][t].push_back(generateNodeSetter<VectorSetter>(r, t, 1, scheme.dirichlet, [] (auto &element, const size_t &n, const size_t &dof, const size_t &s) { return element.temp[n][s]; }));
			}
		}
	}
}

void HeatTransfer::evaluate(SteadyState &scheme, step::Time &time)
{
	setTime(time.current);
	reset(scheme.K, scheme.f, scheme.dirichlet);
	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n", assemble(ActionOperator::Action::ASSEMBLE));
	eslog::info("       = FILL MATRICES                                                  %12.8f s = \n", assemble(ActionOperator::Action::FILL));
	update(scheme.K, scheme.f);
}

void HeatTransfer::dryrun()
{
	if (this->K == nullptr) {
		this->K = new Matrix_Distributed<Matrix_CSR, double>();
		this->K->mapping.elements.resize(info::mesh->elements->eintervals.size());

		this->f = new Vector_Distributed<Vector_Dense, double>();
		this->f->mapping.elements.resize(info::mesh->elements->eintervals.size());
	}
	info::ecf->output.results_selection.flux = !info::ecf->simple_output;
	Assembler::measurements reassemble_time = Assembler::measurements();
	Assembler::measurements   assemble_time = Assembler::measurements();
	Assembler::measurements   solution_time = Assembler::measurements();
	int numreps = 10;
	for(int reps = 0; reps < numreps; reps++) {
		// printf("assemble\n");
		assemble_time += assemble(ActionOperator::Action::ASSEMBLE);
		// printf("reassemble\n");
		reassemble_time += assemble(ActionOperator::Action::REASSEMBLE);
		// printf("solution\n");
		solution_time += assemble(ActionOperator::Action::SOLUTION);
	}

	assemble_time.coreTime /= static_cast<double>(numreps);
	reassemble_time.coreTime /= static_cast<double>(numreps);
	solution_time.coreTime /= static_cast<double>(numreps);

	assemble_time.preprocessTime /= static_cast<double>(numreps);
	reassemble_time.preprocessTime /= static_cast<double>(numreps);
	solution_time.preprocessTime /= static_cast<double>(numreps);

	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n",   assemble_time.preprocessTime);
	std::cout<<"SCALING: "<<assemble_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n",   assemble_time.coreTime);
	std::cout<<"SCALING: "<<assemble_time.coreTime<<std::endl;
	eslog::info("       = SIMD LOOP REASSEMBLE                                           %12.8f s = \n", reassemble_time.preprocessTime);
	std::cout<<"SCALING: "<<reassemble_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP REASSEMBLE                                           %12.8f s = \n", reassemble_time.coreTime);
	std::cout<<"SCALING: "<<reassemble_time.coreTime<<std::endl;
	eslog::info("       = SIMD LOOP SOLUTION                                             %12.8f s = \n", solution_time.preprocessTime);
	std::cout<<"SCALING: "<<solution_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP SOLUTION                                             %12.8f s = \n", solution_time.coreTime);
	std::cout<<"SCALING: "<<solution_time.coreTime<<std::endl;
}

void HeatTransfer::initTemperatureAndThickness()
{
	generateElementExpression<ExternalNodeExpressionWithCoordinates>(etype, elementOps, settings.initial_temperature, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.temp[n][s] = value; });
	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin;
		elementOps[i].push_back(generateElementOperator<InitialTemperatureToNodes>(i, etype[i], procNodes, Results::initialTemperature->data.data()));
		if (info::mesh->dimension == 2) {
			elementOps[i].push_back(generateElementOperator2D<ThicknessToNodes>(i, etype[i], procNodes, Results::thickness->data.data()));
		}
	}
	std::vector<int> bcfilter(info::mesh->boundaryRegions.size(), 0);
	if (settings.init_temp_respect_bc) {
		auto setter = [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.temp[n][s] = value; };

		for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			auto it = configuration.temperature.find(info::mesh->boundaryRegions[r]->name);
			if (it != configuration.temperature.end()) {
				bcfilter[r] = 1;
				for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					auto procNodes = info::mesh->boundaryRegions[r]->nodes->cbegin(t);
					switch (info::mesh->dimension) {
					case 2: boundaryOps[r][t].push_back(generateTypedExpressionNode2D<ExternalNodeExpressionWithCoordinates>(r, t, it->second.evaluator, setter)); break;
					case 3: boundaryOps[r][t].push_back(generateTypedExpressionNode3D<ExternalNodeExpressionWithCoordinates>(r, t, it->second.evaluator, setter)); break;
					}
					boundaryOps[r][t].push_back(generateBoundaryNodeOperator<InitialTemperatureToNodes>(r, t, procNodes, Results::initialTemperature->data.data()));
				}
			}
		}
	}
	assemble(ActionOperator::Action::ASSEMBLE);
	Results::temperature->data = Results::initialTemperature->data;

	dropLastOperators(elementOps);
	dropLastOperators(elementOps);
	if (info::mesh->dimension == 2) {
		dropLastOperators(elementOps);
	}
	if (settings.init_temp_respect_bc) {
		dropLastOperators(bcfilter, boundaryOps);
		dropLastOperators(bcfilter, boundaryOps);
	}
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

template <>
Assembler::measurements HeatTransfer::instantiate2D<HeatTransferElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<HeatTransferDataDescriptor, 1, HeatTransferGPC::POINT1, 2, 0, HeatTransferElementType::NODE>(action, ops, elements);
}

template <>
Assembler::measurements HeatTransfer::instantiate2D<HeatTransferElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<HeatTransferDataDescriptor, 2, HeatTransferGPC::LINE2, 2, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::LINE3, 2, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements HeatTransfer::instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, elements); break;
	default: return measurements();
	};
}

template <>
Assembler::measurements HeatTransfer::instantiate3D<HeatTransferElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<HeatTransferDataDescriptor, 1, HeatTransferGPC::POINT1, 3, 0, HeatTransferElementType::NODE>(action, ops, elements);
}

template <>
Assembler::measurements HeatTransfer::instantiate3D<HeatTransferElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<HeatTransferDataDescriptor, 2, HeatTransferGPC::LINE2, 3, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::LINE3, 3, 1, HeatTransferElementType::EDGE>(action, ops, elements); break;
	default: return measurements();
	}
}

template <>
Assembler::measurements HeatTransfer::instantiate3D<HeatTransferElementType::FACE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 3, 2, HeatTransferElementType::FACE>(action, ops, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements HeatTransfer::instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return loop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return loop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return loop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return loop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return loop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return loop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return loop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return loop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, elements); break;
	default: return measurements();
	};
}

Assembler::measurements HeatTransfer::instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
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
	return measurements();
}

}

