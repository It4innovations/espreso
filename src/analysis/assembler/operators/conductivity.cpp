
#include "conductivity.coordinatesystem.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/heattransfer.generator.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

using namespace espreso;

void HeatTransfer::generateConductivity()
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
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

		elementOps[interval].push_back(generateExpression<ExternalGPsExpression>(interval, etype[interval], mat->density.evaluator,
					[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.density[gp][s] = value; }));
		elementOps[interval].push_back(generateExpression<ExternalGPsExpression>(interval, etype[interval], mat->heat_capacity.evaluator,
					[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatCapacity[gp][s] = value; }));

		const TensorConfiguration &conductivity = mat->thermal_conductivity.values;
		switch (info::mesh->dimension) {
		case 2:
			// [ 0 1 ]
			// [ 2 3 ]
			if (rotate) {
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC: break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][3][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][3][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
									[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][1][s] = element.ecf.conductivity[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][3][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				}

				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator2D<HeatTransferCoordinateSystemCartesian, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator2D<HeatTransferCoordinateSystemCylindric, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				}
				if (settings.reassembling_optimization) {
					addTypedElementStorage2D<HeatTransferElementType::SYMMETRIC_GENERAL>(elementOps[interval], interval, [] (auto &element) { return sizeof(element.cossin); }, [] (auto &element) { return element.cossin; });
				}
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateElementTypedOperator2D<HeatTransferCoordinateSystemApply, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
				elementOps[interval].back()->isconst &= isconst;
			} else {
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][3][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][3][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
									[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][1][s] = element.conductivity[gp][2][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][1][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][2][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression2D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][3][s] = value; }));
					break;
				}
			}
			break;
		case 3:
			// [ 0 1 2 ]
			// [ 3 4 5 ]
			// [ 6 7 8 ]
			if (rotate) {
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC: break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][4][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][8][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][1][s] = element.ecf.conductivity[gp][3][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][2][s] = element.ecf.conductivity[gp][6][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][4][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][5][s] = element.ecf.conductivity[gp][7][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][8][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][3][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][4][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][5][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][6][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][7][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.conductivity[gp][8][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					break;
				}

				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.rotation.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<HeatTransferCoordinateSystemCartesian, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<HeatTransferCoordinateSystemCylindric, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, mat->coordinate_system.center.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementTypedOperator3D<HeatTransferCoordinateSystemSpherical, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
					elementOps[interval].back()->isconst &= isconst;
					break;
				}
				if (settings.reassembling_optimization) {
					addTypedElementStorage3D<HeatTransferElementType::SYMMETRIC_GENERAL>(elementOps[interval], interval, [] (auto &element) { return sizeof(element.cossin); }, [] (auto &element) { return element.cossin; });
				}
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateElementTypedOperator3D<HeatTransferCoordinateSystemApply, HeatTransferElementType::SYMMETRIC_GENERAL>(interval));
				elementOps[interval].back()->isconst &= isconst;
			} else {
				switch (mat->thermal_conductivity.model) {
				case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_ISOTROPIC>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::DIAGONAL:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][4][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][8][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][1][s] = element.conductivity[gp][3][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][2][s] = element.conductivity[gp][6][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][4][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][5][s] = element.conductivity[gp][7][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][8][s] = value; }));
					break;
				case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][0][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][1][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(0, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][2][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][3][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][4][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(1, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][5][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 0).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][6][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 1).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][7][s] = value; }));
					elementOps[interval].push_back(generateTypedExpression3D<ExternalGPsExpression, HeatTransferElementType::SYMMETRIC_GENERAL>(interval, conductivity.get(2, 2).evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.conductivity[gp][8][s] = value; }));
					break;
				}
			}
			break;
		}
	}
}

