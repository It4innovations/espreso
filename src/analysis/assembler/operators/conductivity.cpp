
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

void thermalConductivity(HeatTransfer &module)
{
	bool cooToGP = false;
	for(size_t interval = 0; !cooToGP && interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			cooToGP = true;
		}
	}

	if (module.coords.gp.data == nullptr && cooToGP) {
		module.controller.addInput(module.coords.gp, module.coords.node);
		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
			if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				module.coords.gp.isconst[interval] = true;
				module.coords.gp.update[interval] = -1;
			}
		}
	}
	if (module.coords.gp.data == nullptr) {
		module.controller.prepare(module.coords.gp);
		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				if (info::mesh->dimension == 2) {
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 2, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
				}
				if (info::mesh->dimension == 3) {
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 3, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
				}
			}
		}
	}

	if (cooToGP) {
		if (info::mesh->dimension == 2) {
			module.controller.addInput(module.cooSystem.angle2D, module.coords.gp);
		}
		if (info::mesh->dimension == 3) {
			module.controller.addInput(module.cooSystem.angle3D, module.coords.gp);
		}
	}

	if (info::mesh->dimension == 2) {
		module.controller.addInput(module.cooSystem.angle2D, module.cooSystem.cartesian2D, module.cooSystem.cylindric);
		module.controller.prepare(module.cooSystem.angle2D);
		module.controller.addInput(module.material.conductivity, module.cooSystem.angle2D);
	}
	if (info::mesh->dimension == 3) {
		module.controller.addInput(module.cooSystem.angle3D, module.cooSystem.cartesian3D, module.cooSystem.cylindric, module.cooSystem.spherical);
		module.controller.prepare(module.cooSystem.angle3D);
		module.controller.addInput(module.material.conductivity, module.cooSystem.angle3D);
	}

	module.controller.addInput(module.material.conductivityIsotropic, module.material.model.isotropic);
	module.controller.addInput(module.material.conductivity, module.material.model.diagonal, module.material.model.anisotropic, module.material.model.symmetric2D, module.material.model.symmetric3D);
	module.controller.prepare(module.material.conductivityIsotropic);
	module.controller.prepare(module.material.conductivity);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (info::mesh->dimension == 2) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 1, CopyParameter>(interval, module.controller, module.material.model.isotropic, module.material.conductivityIsotropic)); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopyDiagonal2DConductivity>(interval, module.controller, module.material.model.diagonal, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopySymmetric2DConductivity>(interval, module.controller, module.material.model.symmetric2D, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopyAnisotropic2DConductivity>(interval, module.controller, module.material.model.anisotropic, module.material.conductivity)); break;
			}

			if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					if (module.cooSystem.cartesian2D.externalValues.evaluator[interval]->isset) {
						module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CartesianRotation2D>(interval, module.controller, module.coords.gp, module.cooSystem.cartesian2D, module.cooSystem.angle2D));
						module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotation2D>(interval, module.controller, module.cooSystem.angle2D, module.material.conductivity));
					}
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CylindricalRotation2D>(interval, module.controller, module.coords.gp, module.cooSystem.cylindric, module.cooSystem.angle2D));
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotation2D>(interval, module.controller, module.cooSystem.angle2D, module.material.conductivity));
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					break;
				}
			}
		}
		if (info::mesh->dimension == 3) {
			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, 1, CopyParameter>(interval, module.controller, module.material.model.isotropic, module.material.conductivityIsotropic)); break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopyDiagonal3DConductivity>(interval, module.controller, module.material.model.diagonal, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopySymmetric3DConductivity>(interval, module.controller, module.material.model.symmetric3D, module.material.conductivity)); break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CopyAnisotropic3DConductivity>(interval, module.controller, module.material.model.anisotropic, module.material.conductivity)); break;
			}

			if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					if (module.cooSystem.cartesian3D.externalValues.evaluator[interval]->isset) {
						module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CartesianRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.cartesian3D, module.cooSystem.angle3D));
						module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
					}
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, CylindricalRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.cylindric, module.cooSystem.angle3D));
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, SphericalRotation3D>(interval, module.controller, module.coords.gp, module.cooSystem.spherical, module.cooSystem.angle3D));
					module.elementOps[interval].emplace_back(instantiate<HeatTransfer::NGP, ConductivityRotation3D>(interval, module.controller, module.cooSystem.angle3D, module.material.conductivity));
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

