
#include "heattransfer.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void heatStiffness(HeatTransfer &module)
{
	if (info::mesh->dimension == 2) {
		module.controller.addInput(module.elements.stiffness, module.thickness.gp);
	} else {
		module.thickness.gp.resize();
	}
	module.controller.addInput(module.elements.stiffness, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.material.conductivity, module.material.conductivityIsotropic, module.gradient.xi);
	module.controller.prepare(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(
						instantiate<HeatTransfer::NGP, Stiffness2DHeatIsotropic>(interval, module.controller,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(
						instantiate<HeatTransfer::NGP, Stiffness3DHeatIsotropic>(interval, module.controller,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		} else {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(
						instantiate<HeatTransfer::NGP, Stiffness2DHeat>(interval, module.controller,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(
						instantiate<HeatTransfer::NGP, Stiffness3DHeat>(interval, module.controller,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		}
	}
}

void heatRHS(HeatTransfer &module)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		module.controller.addInput(module.elements.rhs, module.heatSource.gp, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant);
		module.controller.prepare(module.elements.rhs);
		if (module.heatSource.gp.isSet(interval)) {
			module.elementOps[interval].emplace_back(
				instantiate<HeatTransfer::NGP, HeatRHS>(interval, module.controller,
						module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
						module.heatSource.gp,
						module.elements.rhs));
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			bool isSet = module.heatFlow.gp.isSet(r) || module.heatFlux.gp.isSet(r) || module.convection.heatTransferCoeficient.gp.isSet(r);
			if (!isSet) {
				continue;
			}

			module.controller.prepare(module.heatFlow.gp.regions[r], module.heatFlux.gp.regions[r], module.convection.heatTransferCoeficient.gp.regions[r], module.convection.externalTemperature.gp.regions[r]);

			if (info::mesh->dimension == 2) {
				module.controller.addInput(module.elements.boundary.rhs.regions[r], module.thickness.boundary.gp.regions[r]);
			}
			module.controller.addInput(module.q.gp.regions[r], module.heatFlow.gp.regions[r], module.heatFlux.gp.regions[r], module.convection.heatTransferCoeficient.gp.regions[r], module.convection.externalTemperature.gp.regions[r]);
			module.controller.addInput(module.elements.boundary.rhs.regions[r], module.q.gp.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.weight.regions[r]);
			module.controller.prepare(module.q.gp.regions[r], module.elements.boundary.rhs.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryOps[r][interval].emplace_back(instantiate<HeatTransfer::NGP, HeatQ>(r, interval, module.controller,
						info::mesh->boundaryRegions[r]->area,
						module.heatFlow.gp.regions[r],
						module.heatFlux.gp.regions[r],
						module.convection.heatTransferCoeficient.gp.regions[r], module.convection.externalTemperature.gp.regions[r],
						module.q.gp.regions[r]));

				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<HeatTransfer::NGP, HeatRHS2D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.thickness.boundary.gp.regions[r], module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<HeatTransfer::NGP, HeatRHS3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}
