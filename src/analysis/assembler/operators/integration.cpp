
#include "integration.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"


#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void elementIntegration(AX_HeatTransfer &module)
{
	module.integration.jacobiInversion.addInput(module.coords.node);
	module.integration.jacobiInversion.addInput(module.integration.dN);
	module.integration.jacobiInversion.resize();
	module.addParameter(module.integration.jacobiInversion);

	module.integration.jacobiDeterminant.addInput(module.coords.node);
	module.integration.jacobiDeterminant.addInput(module.integration.dN);
	module.integration.jacobiDeterminant.resize();
	module.addParameter(module.integration.jacobiDeterminant);

	module.integration.dND.addInput(module.coords.node);
	module.integration.dND.addInput(module.integration.dN);
	module.integration.dND.resize();
	module.addParameter(module.integration.dND);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ElementJacobian2D>(interval, module.coords.node, module.integration.dN, module.integration.jacobiInversion, module.integration.jacobiDeterminant, module.integration.dND));
		}
		if (info::mesh->dimension == 3) {
			module.actionOps[interval].emplace_back(instantiate<AX_HeatTransfer::NGP, ElementJacobian3D>(interval, module.coords.node, module.integration.dN, module.integration.jacobiInversion, module.integration.jacobiDeterminant, module.integration.dND));
		}
	}
}

}
