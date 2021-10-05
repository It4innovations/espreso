
#include "conductivity.h"
#include "coordinatesystem.h"
#include "copy.h"
#include "gausspoints.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

namespace espreso {

void thermalConductivity(AX_HeatTransfer &module)
{
	if (module.coords.gp.data == nullptr) {
		// we need to compute coordinates in gps
		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				if (info::mesh->dimension == 2) {
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 2, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
					module.controller.addInput(module.cooSystem.angle2D, module.coords.gp);
				}
				if (info::mesh->dimension == 3) {
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 3, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
					module.controller.addInput(module.cooSystem.angle3D, module.coords.gp);
				}
				break;
			}
		}

		// save the memory for cartesian intervals
		module.coords.gp.setConstness(true);
		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				module.coords.gp.isconst[interval] = false;
			}
		}

		module.controller.prepare(module.coords.gp);
	}

	if (info::mesh->dimension == 2) {
		module.controller.addInput(module.cooSystem.angle2D, module.cooSystem.cartesian2D, module.cooSystem.cylindric);
		module.controller.addInput(module.material.conductivity, module.cooSystem.angle2D);
		module.controller.prepare(module.cooSystem.angle2D);
	}
	if (info::mesh->dimension == 3) {
		module.controller.addInput(module.cooSystem.angle3D, module.cooSystem.cartesian3D, module.cooSystem.cylindric, module.cooSystem.spherical);
		module.controller.addInput(module.material.conductivity, module.cooSystem.angle3D);
		module.controller.prepare(module.cooSystem.angle3D);
	}

	module.controller.addInput(module.material.conductivityIsotropic, module.material.model.isotropic);
	module.controller.addInput(module.material.conductivity, module.material.model.diagonal, module.material.model.anisotropic, module.material.model.symmetric2D, module.material.model.symmetric3D);
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

			if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					if (module.cooSystem.cartesian2D.externalValue.evaluator[interval]->isset) {
						module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CartesianRotation2D>(interval, module.controller, module.coords.gp, module.cooSystem.cartesian2D, module.cooSystem.angle2D));
						module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ConductivityRotation2D>(interval, module.controller, module.cooSystem.angle2D, module.material.conductivity));
					}
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CylindricalRotation2D>(interval, module.controller, module.coords.gp, module.cooSystem.cylindric, module.cooSystem.angle2D));
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ConductivityRotation2D>(interval, module.controller, module.cooSystem.angle2D, module.material.conductivity));
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL: break;
				default:
					break;
				}
			}
		}
		if (info::mesh->dimension == 3) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, 1, CopyParameter>(interval, module.controller, module.material.model.isotropic, module.material.conductivityIsotropic)); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyDiagonal3DConductivity>(interval, module.controller, module.material.model.diagonal, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopySymmetric3DConductivity>(interval, module.controller, module.material.model.symmetric2D, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CopyAnisotropic3DConductivity>(interval, module.controller, module.material.model.anisotropic, module.material.conductivity)); break;
			}

			if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					if (module.cooSystem.cartesian3D.externalValue.evaluator[interval]->isset) {
						module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CartesianRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.cartesian3D, module.cooSystem.angle3D));
						module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
					}
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, CylindricalRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.cylindric, module.cooSystem.angle3D));
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL: break;
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, SphericalRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.spherical, module.cooSystem.angle3D));
					module.elementOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
					break;
					break;
				}
			}
		}
		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				module.material.conductivity.update[interval] = -1;
				break;
			default:
				module.material.conductivityIsotropic.update[interval] = -1;
		}
	}
}

}

