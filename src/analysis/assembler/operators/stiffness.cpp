
#include "stiffness.h"

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
				module.actionOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness2DHeatIsotropic>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.actionOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness3DHeatIsotropic>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivityIsotropic,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		} else {
			if (info::mesh->dimension == 2) {
				module.actionOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness2DHeat>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
			if (info::mesh->dimension == 3) {
				module.actionOps[interval].emplace_back(
						instantiate<AX_HeatTransfer::NGP, Stiffness3DHeat>(interval,
								module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.conductivity,
								module.gradient.xi, module.thickness.gp, module.elements.stiffness));
			}
		}
	}
}

}

