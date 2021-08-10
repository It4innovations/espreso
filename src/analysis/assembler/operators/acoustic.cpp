
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
//		module.elements.stiffness.addInput(module.thickness.gp);
	}
	module.elements.stiffness.addInput(module.integration.dND);
	module.elements.stiffness.addInput(module.integration.weight);
	module.elements.stiffness.addInput(module.integration.jacobiDeterminant);
//	module.elements.stiffness.addInput(module.material.conductivity);
//	module.elements.stiffness.addInput(module.gradient.xi);
	module.elements.stiffness.resize();

	module.addParameter(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, Stiffness2DAcoustic>(interval, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.elements.stiffness));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, Stiffness2DAcoustic>(interval, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.elements.stiffness));
		}
	}
}

void acousticMass(AX_Acoustic &module)
{
	if (info::mesh->dimension == 2) {
//		module.elements.stiffness.addInput(module.thickness.gp);
	}
	module.elements.mass.addInput(module.integration.N);
	module.elements.mass.addInput(module.integration.weight);
	module.elements.mass.addInput(module.integration.jacobiDeterminant);
	module.elements.mass.resize();

	module.addParameter(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		module.elementOps[interval].emplace_back(instantiate<AX_Acoustic::NGP, AcousticMass>(interval, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.elements.mass));
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


