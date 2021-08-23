
#include "acoustic.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void acousticStiffness(AX_Acoustic &module)
{
	if (info::mesh->dimension == 2) {
//		module.controller.addInput(module.elements.stiffness, module.thickness.gp);
	}
	module.controller.addInput(module.integration.jacobiInversion, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, Stiffness2DAcoustic>(interval, module.controller, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.elements.stiffness));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, Stiffness2DAcoustic>(interval, module.controller, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.elements.stiffness));
		}
	}
}

void acousticMass(AX_Acoustic &module)
{
	if (info::mesh->dimension == 2) {
//		module.elements.stiffness.addInput(module.thickness.gp);
	}
	module.controller.addInput(module.elements.mass, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticMass>(interval, module.controller, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.elements.mass));
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
			module.controller.addInput(module.q.gp.regions[r], module.normalAcceleration.gp.regions[r]);
			module.controller.addInput(module.elements.boundary.rhs.regions[r], module.q.gp.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.weight.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticQ>(r, interval, module.controller,
						1, // info::mesh->boundaryRegions[r]->area,
						module.normalAcceleration.gp.regions[r],
						module.q.gp.regions[r]));

				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticRHS2D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
//							module.thickness.boundary.gp.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticRHS3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}


