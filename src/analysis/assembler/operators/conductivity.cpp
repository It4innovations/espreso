
#include "conductivity.h"
#include "copy.h"
#include "gausspoints.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

namespace espreso {

void thermalConductivity(HeatTransfer &module)
{
//	bool cooToGP = false;
//	for(size_t interval = 0; !cooToGP && interval < info::mesh->elements->eintervals.size(); ++interval) {
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			cooToGP = true;
//		}
//	}
//
//	if (module.coords.gp.data == nullptr && cooToGP) {
//		module.controller.addInput(module.coords.gp, module.coords.node);
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//			if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//				module.coords.gp.isconst[interval] = true;
//				module.coords.gp.update[interval] = -1;
//			}
//		}
//	}
//	if (module.coords.gp.data == nullptr) {
//		module.controller.prepare(module.coords.gp);
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//				if (info::mesh->dimension == 2) {
//					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 2, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
//				}
//				if (info::mesh->dimension == 3) {
//					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 3, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
//				}
//			}
//		}
//	}
//
//	if (cooToGP) {
//		if (info::mesh->dimension == 2) {
//			module.controller.addInput(module.cooSystem.angle2D, module.coords.gp);
//		}
//		if (info::mesh->dimension == 3) {
//			module.controller.addInput(module.cooSystem.angle3D, module.coords.gp);
//		}
//	}
//
//	if (info::mesh->dimension == 2) {
//		module.controller.addInput(module.cooSystem.angle2D, module.cooSystem.cartesian2D, module.cooSystem.cylindric);
//		module.controller.prepare(module.cooSystem.angle2D);
//		module.controller.addInput(module.material.conductivity, module.cooSystem.angle2D);
//	}
//	if (info::mesh->dimension == 3) {
//		module.controller.addInput(module.cooSystem.angle3D, module.cooSystem.cartesian3D, module.cooSystem.cylindric, module.cooSystem.spherical);
//		module.controller.prepare(module.cooSystem.angle3D);
//		module.controller.addInput(module.material.conductivity, module.cooSystem.angle3D);
//	}
//
//	module.controller.addInput(module.material.conductivityIsotropic, module.material.model.isotropic);
//	module.controller.addInput(module.material.conductivity, module.material.model.diagonal, module.material.model.anisotropic, module.material.model.symmetric2D, module.material.model.symmetric3D);
//	module.controller.prepare(module.material.conductivityIsotropic);
//	module.controller.prepare(module.material.conductivity);
//
//	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//
//		switch (mat->thermal_conductivity.model) {
//		case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityDiagonal, HeatTransferOperator>(interval, module.controller)); break;
//		case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivitySymmetric, HeatTransferOperator>(interval, module.controller)); break;
//		}
//
//		if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
//			switch (mat->coordinate_system.type) {
//			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
//				if (module.cooSystem.cartesian2D.externalValues.evaluator[interval]->isset) {
//					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotationCartesian, HeatTransferOperator>(interval, module.controller));
//					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotationApply, HeatTransferOperator>(interval, module.controller));
//				} else {
//					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotationSkip, HeatTransferOperator>(interval, module.controller));
//				}
//				break;
//			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
//				module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotationCylindric, HeatTransferOperator>(interval, module.controller));
//				module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotationApply, HeatTransferOperator>(interval, module.controller));
//				break;
//			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
//				module.elementOps[interval].emplace_back(instantiate3D<HeatTransfer::NGP, ConductivityRotationSpherical, HeatTransferOperator>(interval, module.controller));
//				module.elementOps[interval].emplace_back(instantiate3D<HeatTransfer::NGP, ConductivityRotationApply, HeatTransferOperator>(interval, module.controller));
//				break;
//			}
//		}
//		switch (mat->thermal_conductivity.model) {
//			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
//				module.material.conductivity.update[interval] = -1;
//				break;
//			default:
//				module.material.conductivityIsotropic.update[interval] = -1;
//		}
//	}
}

}

