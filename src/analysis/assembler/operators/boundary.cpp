
#include "acoustic.h"
#include "heattransfer.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

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

void acousticRHS(AX_Acoustic &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			if (module.normalAcceleration.gp.regions[r].data == NULL) {
				continue;
			}

//			if (info::mesh->dimension == 2) {
//				module.elements.boundary.rhs.regions[r].addInput(module.thickness.boundary.gp.regions[r]);
//			}
			module.q.gp.regions[r].addInput(module.normalAcceleration.gp.regions[r]);
			module.q.gp.regions[r].resize();
			module.addParameter(module.q.gp.regions[r]);

			module.elements.boundary.rhs.regions[r].addInput(module.q.gp.regions[r]);
			module.elements.boundary.rhs.regions[r].addInput(module.integration.boundary.jacobian.regions[r]);
			module.elements.boundary.rhs.regions[r].addInput(module.integration.boundary.weight.regions[r]);
			module.elements.boundary.rhs.regions[r].resize();
			module.addParameter(module.elements.boundary.rhs.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticQ>(r, interval,
						1, // info::mesh->boundaryRegions[r]->area,
						module.normalAcceleration.gp.regions[r],
						module.q.gp.regions[r]));

				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticRHS2D>(r, interval,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
//							module.thickness.boundary.gp.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticRHS3D>(r, interval,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}
