
#include "flux.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"
#include "config/ecf/material/material.h"

namespace espreso {

void outputFlux(AX_HeatTransfer &module)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
		if (mat->thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
			if (info::mesh->dimension == 2) {
				module.elementRes[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, OutputFluxIsotropic2D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivityIsotropic, module.flux.output));
			}
			if (info::mesh->dimension == 3) {
				module.elementRes[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, OutputFluxIsotropic3D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivityIsotropic, module.flux.output));
			}
		} else {
			if (info::mesh->dimension == 2) {
				module.elementRes[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, OutputFlux2D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivity, module.flux.output));
			}
			if (info::mesh->dimension == 3) {
				module.elementRes[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, OutputFlux3D>(interval, module.controller, module.integration.dND, module.temp.node, module.material.conductivity, module.flux.output));
			}
		}
	}
}

}
