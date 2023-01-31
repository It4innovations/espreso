
#include "heat.h"
#include "heattransfer.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/coordinatesystem.h"
#include "analysis/assembler/operators/heattransfer.stiffness.h"

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

NodeData* Heat::Results::temperature = nullptr;
NodeData* Heat::Results::initialTemperature = nullptr;
ElementData* Heat::Results::translationMotion = nullptr;
ElementData* Heat::Results::gradient = nullptr;
ElementData* Heat::Results::flux = nullptr;

Heat::Heat(Heat *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{

}

void Heat::initParameters()
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

bool Heat::initTemperature()
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
//	// TODO: fix
////	fromExpression(*this, temp.initial.node, temp.initial.node.externalValues);
//	_evaluate();
//
////	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
////		Heat::insertParameters(it->second.evaluator);
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

template<int etype, int i, int j, int index>
void pushConductivity2D(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, const MaterialConfiguration *mat, bool rotate)
{
	if (rotate) {
		ops[interval].push_back(generateTypedExpression2D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][index][s] = value; }));
	} else {
		ops[interval].push_back(generateTypedExpression2D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][index][s] = value; }));
	}
}

template<int etype, int i, int j, int index0, int index1>
void pushConductivity2D2(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, const MaterialConfiguration *mat, bool rotate)
{
	if (rotate) {
		ops[interval].push_back(generateTypedExpression2D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][index0][s] = element.ecf.conductivity[gp][index1][s] = value; }));
	} else {
		ops[interval].push_back(generateTypedExpression2D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][index0][s] = element.conductivity[gp][index1][s] = value; }));
	}
}

template<int i, int j, int index>
void pushConductivity2D(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, int etype, const MaterialConfiguration *mat, bool rotate)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL: pushConductivity2D<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index>(ops, interval, mat, rotate); break;
//	case HeatTransferElementType::ASYMMETRIC_GENERAL: pushConductivity2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate);
	}
}

template<int i, int j, int index0, int index1>
void pushConductivity2D2(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, int etype, const MaterialConfiguration *mat, bool rotate)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL: pushConductivity2D2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate); break;
//	case HeatTransferElementType::ASYMMETRIC_GENERAL: pushConductivity2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate);
	}
}

template<int etype, int i, int j, int index>
void pushConductivity3D(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, const MaterialConfiguration *mat, bool rotate)
{
	if (rotate) {
		ops[interval].push_back(generateTypedExpression3D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][index][s] = value; }));
	} else {
		ops[interval].push_back(generateTypedExpression3D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][index][s] = value; }));
	}
}

template<int etype, int i, int j, int index0, int index1>
void pushConductivity3D2(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, const MaterialConfiguration *mat, bool rotate)
{
	if (rotate) {
		ops[interval].push_back(generateTypedExpression3D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][index0][s] = element.ecf.conductivity[gp][index1][s] = value; }));
	} else {
		ops[interval].push_back(generateTypedExpression3D<ExternalExpression, etype>(interval, mat->thermal_conductivity.values.get(i, j).evaluator,
				[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][index0][s] = element.conductivity[gp][index1][s] = value; }));
	}
}

template<int i, int j, int index>
void pushConductivity3D(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, int etype, const MaterialConfiguration *mat, bool rotate)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL: pushConductivity3D<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index>(ops, interval, mat, rotate); break;
//	case HeatTransferElementType::ASYMMETRIC_GENERAL: pushConductivity2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate);
	}
}

template<int i, int j, int index0, int index1>
void pushConductivity3D2(std::vector<std::vector<ActionOperator*> > &ops, size_t interval, int etype, const MaterialConfiguration *mat, bool rotate)
{
	switch (etype) {
	case HeatTransferElementType::SYMMETRIC_GENERAL: pushConductivity3D2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate); break;
//	case HeatTransferElementType::ASYMMETRIC_GENERAL: pushConductivity2<HeatTransferElementType::SYMMETRIC_GENERAL, i, j, index0, index1>(ops, interval, mat, rotate);
	}
}

void Heat::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");
//	initNames();

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);
	etype.resize(info::mesh->elements->eintervals.size(), 0);
	std::vector<int> cooToGP(info::mesh->elements->eintervals.size(), 0);

	initParameters();

	eslog::info(" ============================================================================================= \n");
	bool correct = true;
//	correct &= initTemperature();

	if (step::step.loadstep == 0) {
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
					correct &= checkMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x);
					correct &= checkMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y);
					correct &= checkMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					correct &= checkMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x);
					correct &= checkMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y);
					correct &= checkMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x);
					correct &= checkMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x);
					correct &= checkMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			correct &= checkMaterialParameter(mat->name, "DENSITY", mat->density);
			correct &= checkMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2));
					correct &= checkMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0));
					correct &= checkMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1));
					correct &= checkMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0));
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:   etype[interval] = HeatTransferElementType::SYMMETRIC_ISOTROPIC; break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:    etype[interval] = HeatTransferElementType::SYMMETRIC_GENERAL  ; break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   etype[interval] = HeatTransferElementType::SYMMETRIC_GENERAL  ; break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: etype[interval] = HeatTransferElementType::ASYMMETRIC_GENERAL ; break;
			}
			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				cooToGP[interval] |= mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
			}
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);
		eslog::info(" ============================================================================================= \n");
	}

	// we cannot generate functions since 'etype' is filled
	generateBaseFunctions(etype, elementOps);
//	generateElementOperators<CoordinatesToElementNodes2>(etype, elementOps, info::mesh->elements->nodes->cbegin());

//		if (info::mesh->dimension == 2) {
//			correct &= examineElementParameter("THICKNESS", settings.thickness, thickness.gp.externalValues);
//			fromExpression2D(*this, thickness.gp, thickness.gp.externalValues, [] (auto &element) { return element.ecf.thickness; });
//		}

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (cooToGP[interval]) {
			elementOps[interval].push_back(generateElementOperator<CoordinatesToElementNodesAndGPs>(interval, etype[interval], info::mesh->elements->nodes->cbegin()));
		} else {
			elementOps[interval].push_back(generateElementOperator<CoordinatesToElementNodes2>(interval, etype[interval], info::mesh->elements->nodes->cbegin()));
		}

		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];

		bool isconst = true, rotate = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
		if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				if (info::mesh->dimension == 2) {
					rotate &= mat->coordinate_system.rotation.z.isset;
				}
				if (info::mesh->dimension == 3) {
					rotate &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
				}
			}
		}

		switch (info::mesh->dimension) {
		case 2:
			// [ 0 1 ]
			// [ 2 3 ]
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				elementOps[interval].push_back(generateTypedExpression2D<ExternalExpression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, mat->thermal_conductivity.values.get(0, 0).evaluator,
						[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][s] = value; }));
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				pushConductivity2D<0, 0, 0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D<1, 1, 3>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				pushConductivity2D <0, 0,    0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D2<0, 1, 1, 2>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D <1, 1,    3>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				pushConductivity2D<0, 0, 0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D<0, 1, 1>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D<1, 0, 2>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity2D<1, 1, 3>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			}

			if (rotate) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator2D<CoordinateSystemCartesian, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator2D<CoordinateSystemCylindric, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				}
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateElementTypedOperator2D<CoordinateSystemApply, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
				elementOps[interval].back()->isconst &= isconst;
			}
			break;
		case 3:
			// [ 0 1 2 ]
			// [ 3 4 5 ]
			// [ 6 7 8 ]
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, mat->thermal_conductivity.values.get(0, 0).evaluator,
						[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][s] = value; }));
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				pushConductivity3D<0, 0, 0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<1, 1, 4>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<2, 2, 8>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				pushConductivity3D <0, 0,    0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D2<0, 1, 1, 3>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D2<0, 2, 2, 6>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D <1, 1,    4>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D2<1, 2, 5, 7>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D <2, 2,    8>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				pushConductivity3D<0, 0, 0>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<0, 1, 1>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<0, 2, 2>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<1, 0, 3>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<1, 1, 4>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<1, 2, 5>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<2, 0, 6>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<2, 1, 7>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				pushConductivity3D<2, 2, 8>(elementOps, interval, etype[interval], mat, rotate); isconst &= elementOps[interval].back()->isconst;
				break;
			}
			if (rotate) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<CoordinateSystemCartesian, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<CoordinateSystemCylindric, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angle[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<CoordinateSystemSpherical, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				}
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateElementTypedOperator3D<CoordinateSystemApply, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
				elementOps[interval].back()->isconst &= isconst;
			}
			break;
		}
	}

	generateElementOperators<Integration>(etype, elementOps);
//	volume();

//	gradient.xi.resize(1);
//	controller.prepare(gradient.xi);
	elements.stiffness.setConstness(false);
	elements.stiffness.resize();
	generateElementOperators<HeatTransferStiffness2>(etype, elementOps, elements.stiffness);

//	if (configuration.heat_source.size()) {
//		correct &= examineElementParameter("HEAT SOURCE", configuration.heat_source, heatSource.gp.externalValues);
//		fromExpression(*this, heatSource.gp, heatSource.gp.externalValues);
//	}
//	if (configuration.heat_flow.size()) {
//		correct &= examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, heatFlow.gp.externalValues);
//		fromExpression(*this, heatFlow.gp, heatFlow.gp.externalValues);
//	}
//	if (configuration.heat_flux.size()) {
//		correct &= examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, heatFlux.gp.externalValues);
//		fromExpression(*this, heatFlux.gp, heatFlux.gp.externalValues);
//	}
//	RHS(*this);
//
//	outputGradient(*this);
//	outputFlux(*this);
//
//	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
//	eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
//	eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize<Heat, HeatOperator>());
//	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
//	if (correct) {
//		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
//	} else {
//		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
//	}
//	eslog::info(" ============================================================================================= \n");
}

void Heat::connect(SteadyState &scheme)
{
//	addFiller(*this, scheme);
}

void Heat::evaluate(SteadyState &scheme)
{
//	controller.setUpdate();
//	reset(scheme.K, scheme.f, scheme.dirichlet);
////	iterate();
////	printVersions();
	eslog::info("       = SIMD LOOP                                                      %12.8f s = \n", assemble());
	eslog::info("       = SIMD LOOP                                                      %12.8f s = \n", assemble());
	eslog::info("       = SIMD LOOP                                                      %12.8f s = \n", assemble());
////	fill();
//	update(scheme.K, scheme.f);
//	controller.resetUpdate();
}

void Heat::volume()
{
	std::vector<double> volume(info::mesh->elements->eintervals.size());

	generateElementOperators<Volume>(etype, elementOps, volume);
	assemble();
	printElementVolume(volume);
//	printBoundarySurface();
	dropLastOperators(elementOps);
}

void Heat::updateSolution(SteadyState &scheme)
{
//	scheme.x->storeTo(Results::temperature->data);
//	results(); // do we need an update mechanism?
//	temp.node.setUpdate(1);
}

template <int etype>
double Heat::typedInstance(size_t interval, int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TRIANGLE3): return simdloop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::TRIANGLE6): return simdloop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE4):   return simdloop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::SQUARE8):   return simdloop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(ops, elements); break;
		}
		break;
	case 3:
		switch (code) {
		case static_cast<size_t>(Element::CODE::TETRA4):    return loop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::TETRA10):   return loop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID5):  return loop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::PYRAMID13): return loop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA6):   return loop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::PRISMA15):  return loop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA8):     return loop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(ops, elements); break;
		case static_cast<size_t>(Element::CODE::HEXA20):    return loop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(ops, elements); break;
		}
		break;
	}
	return 0;
}

double Heat::instantiate(size_t interval, int code, const std::vector<ActionOperator*> &ops, esint elements)
{
	switch (etype[interval]) {
	case HeatTransferElementType::SYMMETRIC_ISOTROPIC: return typedInstance<HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, code, ops, elements);
	case HeatTransferElementType::SYMMETRIC_GENERAL:   return typedInstance<HeatTransferElementType::SYMMETRIC_GENERAL  >(interval, code, ops, elements);
	}
	return 0;
}

void Heat::initNames()
{

}

void Heat::printVersions()
{

}
