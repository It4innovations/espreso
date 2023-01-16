
#include "acoustic.forces.h"
#include "acoustic.stiffness.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "analysis/assembler/operators/boundarynormals.h"

#include <iostream>

namespace espreso {

void stiffness(Acoustic &module)
{
	if (info::mesh->dimension == 2) {
//		module.controller.addInput(module.elements.stiffness, module.thickness.gp);
	}
	module.controller.addInput(module.elements.stiffness, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant);
	module.controller.prepare(module.elements.stiffness);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, Stiffness2DAcoustic>(interval, module.controller, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.elements.stiffness));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, Stiffness3DAcoustic>(interval, module.controller, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.elements.stiffness));
		}
	}
}

void mass(Acoustic &module)
{
	if (info::mesh->dimension == 2) {
//		module.elements.stiffness.addInput(module.thickness.gp);
	}
	module.controller.addInput(module.elements.mass, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.material.speed_of_sound);
	module.controller.prepare(module.elements.mass);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, AcousticMass>(interval, module.controller, module.integration.N, module.integration.weight, module.integration.jacobiDeterminant, module.material.density, module.material.speed_of_sound ,module.elements.mass));
	}
}

void boundaryMass(Acoustic &module)
{
}

void RHS(Acoustic &module)
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (module.monopoleSource.gp.isSet(interval)) {
			module.controller.addInput(module.elements.monopole, module.monopoleSource.gp,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.material.density);
			module.controller.prepare(module.elements.monopole);

			module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, AcousticRHS2D>(interval, module.controller,
								module.integration.N, module.integration.weight, module.integration.jacobiDeterminant,
								module.monopoleSource.gp,
								module.elements.monopole));
		}

		if (module.dipoleSource.gp.isSet(interval)) {
			module.controller.addInput(module.elements.dipole, module.dipoleSource.gp, module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant);
			module.controller.prepare(module.elements.dipole);

			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, AcousticDipole2D>(interval, module.controller,
									module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density,
									module.dipoleSource.gp,
									module.elements.dipole));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(instantiate<Acoustic::NGP, AcousticDipole3D>(interval, module.controller,
									module.integration.dND, module.integration.weight, module.integration.jacobiDeterminant,
									module.material.density,
									module.dipoleSource.gp,
									module.elements.dipole));
			}
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			if (!module.normalAcceleration.gp.isSet(r) && !module.impedance.gp.isSet(r) && !module.acceleration.gp.isSet(r)) {
				continue;
			}

			bool impedance = module.impedance.gp.isSet(r);
			module.controller.prepare(module.normalAcceleration.gp.regions[r], module.impedance.gp.regions[r]);

			module.controller.addInput(module.q.gp.regions[r], module.normalAcceleration.gp.regions[r], module.impedance.gp.regions[r]);
			module.controller.addInput(module.elements.boundary.rhs.regions[r], module.q.gp.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.weight.regions[r]);
			module.controller.prepare(module.q.gp.regions[r], module.elements.boundary.rhs.regions[r]);
			if (impedance) {
				module.controller.addInput(module.elements.boundary.mass.regions[r], module.impedance.gp.regions[r], module.integration.boundary.N.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.weight.regions[r]);
				module.controller.prepare(module.elements.boundary.mass.regions[r]);
			}

			module.controller.addInput(module.normals.gp.regions[r], module.integration.boundary.dN.regions[r], module.coords.boundary.node.regions[r]);
			module.controller.prepare(module.normals.gp.regions[r]);

			module.controller.addInput(module.proj_acceleration.gp.regions[r], module.acceleration.gp.regions[r], module.normals.gp.regions[r], module.integration.boundary.N.regions[r], module.integration.boundary.jacobian.regions[r], module.integration.boundary.weight.regions[r]);
			module.controller.prepare(module.proj_acceleration.gp.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticQ>(r, interval, module.controller,
						1, // info::mesh->boundaryRegions[r]->area,
						module.normalAcceleration.gp.regions[r],
						module.impedance.gp.regions[r],
						module.q.gp.regions[r]));

				if (impedance) {
					module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticsBoundaryMass>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.impedance.gp.regions[r],
							module.elements.boundary.mass.regions[r]));
				}

				if (info::mesh->dimension == 2) {
					if (module.acceleration.gp.isSet(r)) {
						module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, BoundaryNormal2D>(r, interval, module.controller,
								module.integration.boundary.dN.regions[r],
								module.coords.boundary.node.regions[r],
								module.normals.gp.regions[r]));
						module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticAcceleration2D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.normals.gp.regions[r],
							module.acceleration.gp.regions[r],
							module.proj_acceleration.gp.regions[r]
						));
					}
					module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticRHS2D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					if (module.acceleration.gp.isSet(r)) {
						module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, BoundaryNormal3D>(r, interval, module.controller,
								module.integration.boundary.dN.regions[r],
								module.coords.boundary.node.regions[r],
								module.normals.gp.regions[r]));
						module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticAcceleration3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.normals.gp.regions[r],
							module.acceleration.gp.regions[r],
							module.proj_acceleration.gp.regions[r]
						));
					}

					module.boundaryOps[r][interval].emplace_back(instantiate<Acoustic::NGP, AcousticRHS3D>(r, interval, module.controller,
							module.integration.boundary.N.regions[r], module.integration.boundary.weight.regions[r], module.integration.boundary.jacobian.regions[r],
							module.q.gp.regions[r],
							module.elements.boundary.rhs.regions[r]));
				}
			}
		}
	}
}

}


