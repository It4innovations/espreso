
#include "flux.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"
#include "config/ecf/material/material.h"

namespace espreso {

void outputFlux(HeatTransfer &module)
{
//	if (HeatTransfer::Results::flux) {
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//			if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
//				if (info::mesh->dimension == 2) {
//					module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputFluxIsotropic2D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivityIsotropic, HeatTransfer::Results::flux));
//				}
//				if (info::mesh->dimension == 3) {
//					module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputFluxIsotropic3D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivityIsotropic, HeatTransfer::Results::flux));
//				}
//			} else {
//				if (info::mesh->dimension == 2) {
//					module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputFlux2D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivity, HeatTransfer::Results::flux));
//				}
//				if (info::mesh->dimension == 3) {
//					module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputFlux3D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivity, HeatTransfer::Results::flux));
//				}
//			}
//		}
//	}
}

}
