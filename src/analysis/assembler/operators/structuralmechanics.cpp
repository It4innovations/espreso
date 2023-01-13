
#include "structuralmechanics.h"
#include "analysis/assembler/module/structuralmechanics.h"

#include "analysis/assembler/operator.hpp"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

void stiffness(StructuralMechanics &module)
{
	if (info::mesh->dimension == 2) {
		module.controller.addInput(module.elements.stiffness, module.thickness.gp);
	} else {
		module.thickness.gp.resize();
	}
	module.controller.addInput(module.elements.stiffness, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.material.elasticity2D, module.material.elasticity2DAxisymm, module.material.elasticity3D);
	module.controller.prepare(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			switch (module.settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
				module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, Stiffness2DPlane>(interval, module.controller,
									module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.elasticity2D, module.thickness.gp, module.elements.stiffness));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, Stiffness2DPlaneWithThickness>(interval, module.controller,
									module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.elasticity2D, module.thickness.gp, module.elements.stiffness));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				break;
			}
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, Stiffness3DElasticity>(interval, module.controller,
									module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.elasticity3D, module.thickness.gp, module.elements.stiffness));
			module.controller.prepare(module.material.elasticity3D);
		}
	}
}

void RHS(StructuralMechanics &module)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		module.controller.addInput(module.elements.rhs, module.acceleration.gp, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.thickness.gp);
		module.controller.prepare(module.elements.rhs);
		if (module.acceleration.gp.isSet(interval)) {
			switch (info::mesh->dimension) {
			case 2:
				module.elementOps[interval].emplace_back(
						instantiate<StructuralMechanics::NGP, Acceleration2D>(interval, module.controller,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.density, module.thickness.gp,
								module.acceleration.gp,
								module.elements.rhs));
				break;
			case 3:
				module.elementOps[interval].emplace_back(
						instantiate<StructuralMechanics::NGP, Acceleration3D>(interval, module.controller,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.density, module.thickness.gp,
								module.acceleration.gp,
								module.elements.rhs));
				break;
			}

		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			bool isSet = module.normalPressure.gp.isSet(r);
			if (!isSet) {
				continue;
			}

			module.controller.prepare(module.normalPressure.gp.regions[r]);

			if (info::mesh->dimension == 2) {
				module.controller.addInput(module.elements.boundary.rhs.regions[r], module.thickness.boundary.gp.regions[r]);
			}
			module.controller.addInput(module.elements.boundary.rhs.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.normalPressure.gp.regions[r]);
			module.controller.prepare(module.elements.boundary.rhs.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
//				module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, HeatQ>(r, interval, module.controller,
//						info::mesh->boundaryRegions[r]->area,
//						module.heatFlow.gp.regions[r],
//						module.heatFlux.gp.regions[r],
//						module.convection.heatTransferCoeficient.gp.regions[r], module.convection.externalTemperature.gp.regions[r],
//						module.q.gp.regions[r]));
//
				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressure2D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
							module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r],
							module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressure3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
							module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r],
							module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}

