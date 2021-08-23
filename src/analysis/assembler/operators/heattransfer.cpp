
#include "heattransfer.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void heatStiffness(AX_HeatTransfer &module)
{
	if (info::mesh->dimension == 2) {
		module.elements.stiffness.addInput(module.thickness.gp);
	} else {
		module.thickness.gp.resize();
	}
	module.elements.stiffness.addInput(module.integration.dND);
	module.elements.stiffness.addInput(module.integration.weight);
	module.elements.stiffness.addInput(module.integration.jacobiDeterminant);
	module.elements.stiffness.addInput(module.material.conductivity);
	module.elements.stiffness.addInput(module.gradient.xi);
	module.elements.stiffness.resize();

	module.addParameter(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness2DHeatIsotropic>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness3DHeatIsotropic>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		} else {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness2DHeat>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness3DHeat>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		}
	}
}

void heatRHS(AX_HeatTransfer &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			if (
					module.heatFlow.gp.regions[r].data == NULL &&
					module.heatFlux.gp.regions[r].data == NULL &&
					module.convection.heatTransferCoeficient.gp.regions[r].data == NULL
					) {
				continue;
			}

			if (module.heatFlow.gp.regions[r].data == NULL) {
				module.heatFlow.gp.regions[r].resize();
				module.addParameter(module.heatFlow.gp.regions[r]);
			}
			if (module.heatFlux.gp.regions[r].data == NULL) {
				module.heatFlux.gp.regions[r].resize();
				module.addParameter(module.heatFlux.gp.regions[r]);
			}
			if (module.convection.heatTransferCoeficient.gp.regions[r].data == NULL) {
				module.convection.heatTransferCoeficient.gp.regions[r].resize();
				module.convection.externalTemperature.gp.regions[r].resize();
				module.addParameter(module.convection.heatTransferCoeficient.gp.regions[r]);
				module.addParameter(module.convection.externalTemperature.gp.regions[r]);
			}

			if (info::mesh->dimension == 2) {
				module.elements.boundary.rhs.regions[r].addInput(module.thickness.boundary.gp.regions[r]);
			}
			module.q.gp.regions[r].addInput(module.heatFlow.gp.regions[r]);
			module.q.gp.regions[r].addInput(module.heatFlux.gp.regions[r]);
			module.q.gp.regions[r].addInput(module.convection.heatTransferCoeficient.gp.regions[r]);
			module.q.gp.regions[r].addInput(module.convection.externalTemperature.gp.regions[r]);
			module.q.gp.regions[r].resize();
			module.addParameter(module.q.gp.regions[r]);

			module.elements.boundary.rhs.regions[r].addInput(module.q.gp.regions[r]);
			module.elements.boundary.rhs.regions[r].addInput(module.integration.boundary.jacobian.regions[r]);
			module.elements.boundary.rhs.regions[r].addInput(module.integration.boundary.weight.regions[r]);
			module.elements.boundary.rhs.regions[r].resize();
			module.addParameter(module.elements.boundary.rhs.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryOps[r][interval].emplace_back(instantiate<AX_HeatTransfer::NGP, HeatQ>(r, interval,
						info::mesh->boundaryRegions[r]->area,
						module.heatFlow.gp.regions[r],
						module.heatFlux.gp.regions[r],
						module.convection.heatTransferCoeficient.gp.regions[r], module.convection.externalTemperature.gp.regions[r],
						module.q.gp.regions[r]));

				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_HeatTransfer::NGP, HeatRHS2D>(r, interval,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.thickness.boundary.gp.regions[r], module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_HeatTransfer::NGP, HeatRHS3D>(r, interval,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}
