
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
		if (module.settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC) {
			module.controller.addInput(module.elements.stiffness, module.coords.gp);
		}
	} else {
		module.thickness.gp.resize();
	}
	module.controller.addInput(module.elements.stiffness, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.material.elasticityPlane, module.material.elasticityAxisymm, module.material.elasticity3D);
	module.controller.prepare(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			switch (module.settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
				module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, StiffnessPlane>(interval, module.controller,
									module.integration.N, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.coords.gp,
									module.material.elasticityPlane, module.thickness.gp, module.elements.stiffness));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, StiffnessPlaneWithThickness>(interval, module.controller,
									module.integration.N, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.coords.gp,
									module.material.elasticityPlane, module.thickness.gp, module.elements.stiffness));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, StiffnessAxisymmetric>(interval, module.controller,
									module.integration.N, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.coords.gp,
									module.material.elasticityAxisymm, module.thickness.gp, module.elements.stiffness));
				break;
			}
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(
									instantiate<StructuralMechanics::NGP, Stiffness3DElasticity>(interval, module.controller,
									module.integration.N, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.coords.gp,
									module.material.elasticity3D, module.thickness.gp, module.elements.stiffness));
			module.controller.prepare(module.material.elasticity3D);
		}
	}
}

void RHS(StructuralMechanics &module)
{
	if (module.acceleration.gp.isSet() || module.angularVevocity.gp.isSet()) {
		module.controller.addInput(module.elements.rhs, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.thickness.gp);
	}
	if (module.acceleration.gp.isSet()) {
		module.controller.addInput(module.elements.rhs, module.acceleration.gp);
	}
	if (module.angularVevocity.gp.isSet()) {
		module.controller.addInput(module.elements.rhs, module.angularVevocity.gp, module.coords.gp);
	}
	module.controller.prepare(module.elements.rhs);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.acceleration.gp.isSet(interval)) {
			switch (info::mesh->dimension) {
			case 2:
				switch (module.settings.element_behaviour) {
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AccelerationPlane>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.coords.gp, module.thickness.gp,
									module.acceleration.gp,
									module.elements.rhs));
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AccelerationPlaneWithThickness>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.coords.gp, module.thickness.gp,
									module.acceleration.gp,
									module.elements.rhs));
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AccelerationAxisymmetric>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.coords.gp, module.thickness.gp,
									module.acceleration.gp,
									module.elements.rhs));
					break;
				}
				break;
			case 3:
				module.elementOps[interval].emplace_back(
						instantiate<StructuralMechanics::NGP, Acceleration3D>(interval, module.controller,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.density, module.coords.gp, module.thickness.gp,
								module.acceleration.gp,
								module.elements.rhs));
				break;
			}
		}
		if (module.angularVevocity.gp.isSet(interval)) {
			switch (info::mesh->dimension) {
			case 2:
				switch (module.settings.element_behaviour) {
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AngularVelocityPlane>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.thickness.gp,
									module.coords.gp, module.angularVevocity.gp,
									module.elements.rhs));
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AngularVelocityPlaneWithThickness>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.thickness.gp,
									module.coords.gp, module.angularVevocity.gp,
									module.elements.rhs));
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
					module.elementOps[interval].emplace_back(
							instantiate<StructuralMechanics::NGP, AngularVelocityAxisymmetric>(interval, module.controller,
									module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density, module.thickness.gp,
									module.coords.gp, module.angularVevocity.gp,
									module.elements.rhs));
					break;
				}
				break;
			case 3:
				module.elementOps[interval].emplace_back(
						instantiate<StructuralMechanics::NGP, AngularVelocity3D>(interval, module.controller,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.density, module.thickness.gp,
								module.coords.gp, module.angularVevocity.gp,
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
				if (module.settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC) {
					module.controller.addInput(module.elements.boundary.rhs.regions[r], module.coords.boundary.gp.regions[r]);
				}
			}
			module.controller.addInput(module.elements.boundary.rhs.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.normalPressure.gp.regions[r]);
			module.controller.prepare(module.elements.boundary.rhs.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				if (info::mesh->dimension == 2) {
					switch (module.settings.element_behaviour) {
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
						module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressurePlane>(r, interval, module.controller,
										module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
										module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.coords.boundary.gp.regions[r],
										module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
										module.elements.boundary.rhs.regions[r]));
						break;
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
						module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressurePlaneWithThickness>(r, interval, module.controller,
										module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
										module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.coords.boundary.gp.regions[r],
										module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
										module.elements.boundary.rhs.regions[r]));
						break;
					case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
						module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressureAxisymmetric>(r, interval, module.controller,
										module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
										module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.coords.boundary.gp.regions[r],
										module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
										module.elements.boundary.rhs.regions[r]));
						break;
					}
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<StructuralMechanics::NGP, NormalPressure3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r],
							module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r], module.coords.boundary.gp.regions[r],
							module.thickness.boundary.gp.regions[r], module.normalPressure.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}

