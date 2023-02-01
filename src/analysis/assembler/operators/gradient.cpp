
#include "gradient.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

void outputGradient(HeatTransfer &module)
{
//	if (HeatTransfer::Results::gradient) {
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			if (info::mesh->dimension == 2) {
//				module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputGradient2D>(interval, module.controller, module.integration.dND, module.temp.node, HeatTransfer::Results::gradient));
//			}
//			if (info::mesh->dimension == 3) {
//				module.elementRes[interval].emplace_back(instantiate<HeatTransfer::NGP, OutputGradient3D>(interval, module.controller, module.integration.dND, module.temp.node, HeatTransfer::Results::gradient));
//			}
//		}
//	}
}

}
