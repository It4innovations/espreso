
#include "conductivity.h"
#include "copy.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

namespace espreso {

void thermalConductivity(AX_HeatTransfer &module)
{
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN: break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL: module.controller.addInput(i, module.material.conductivity, info::mesh->nodes->coordinates); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: module.controller.addInput(i, module.material.conductivity, info::mesh->nodes->coordinates); break;
			break;
		}
	}

	module.controller.addInput(module.material.conductivityIsotropic, module.material.model.isotropic);
	module.controller.addInput(module.material.conductivity, module.material.model.diagonal, module.material.model.anisotropic, module.cooSystem.spherical, module.cooSystem.cylindric);
	if (info::mesh->dimension == 2) {
		module.controller.addInput(module.material.conductivity, module.material.model.symmetric2D, module.cooSystem.cartesian2D);
	}
	if (info::mesh->dimension == 3) {
		module.controller.addInput(module.material.conductivity, module.material.model.symmetric3D, module.cooSystem.cartesian3D);
	}

	module.controller.prepare(module.material.conductivityIsotropic);
	module.controller.prepare(module.material.conductivity);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (info::mesh->dimension == 2) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, CopyParameter>(interval, module.controller, module.material.model.isotropic, module.material.conductivityIsotropic)); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyDiagonal2DConductivity>(interval, module.controller, module.material.model.diagonal, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopySymmetric2DConductivity>(interval, module.controller, module.material.model.symmetric2D, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyAnisotropic2DConductivity>(interval, module.controller, module.material.model.anisotropic, module.material.conductivity)); break;
			}

//			switch (mat->coordinate_system.type) {
//			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<typename Assembler::NGP>(Cartesian2DCoordinateSystem(module.coords.gp, module.cooSystem.cartesian2D, module.material.conductivity, interval)); break;
//			case CoordinateSystemConfiguration::TYPE::SPHERICAL: break;
//			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<typename Assembler::NGP>(Cylindrical2DCoordinateSystem(module.coords.gp, module.cooSystem.cylindric, module.material.conductivity, interval)); break;
//				break;
//			}
		}
		if (info::mesh->dimension == 3) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, CopyParameter>(interval, module.controller, module.material.model.isotropic, module.material.conductivityIsotropic)); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyDiagonal3DConductivity>(interval, module.controller, module.material.model.diagonal, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopySymmetric3DConductivity>(interval, module.controller, module.material.model.symmetric2D, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyAnisotropic3DConductivity>(interval, module.controller, module.material.model.anisotropic, module.material.conductivity)); break;
			}

//			switch (mat->coordinate_system.type) {
//			case CoordinateSystemConfiguration::TYPE::CARTESIAN: iterate_elements_gps<typename Assembler::NGP>(Cartesian3DCoordinateSystem(module.coords.gp, module.cooSystem.cartesian3D, module.material.conductivity, interval)); break;
//			case CoordinateSystemConfiguration::TYPE::SPHERICAL: iterate_elements_gps<typename Assembler::NGP>(Spherical3DCoordinateSystem(module.coords.gp, module.cooSystem.spherical, module.material.conductivity, interval)); break;
//			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: iterate_elements_gps<typename Assembler::NGP>(Cylindrical3DCoordinateSystem(module.coords.gp, module.cooSystem.cylindric, module.material.conductivity, interval)); break;
//				break;
//			}
		}
		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.material.conductivity.update[interval] = -1; break;
			default: module.material.conductivityIsotropic.update[interval] = -1; break;
		}
	}
}

}

