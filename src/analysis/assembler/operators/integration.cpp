
#include "integration.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template <class Module>
void _elementIntegration(Module &module)
{
	module.controller.addInput(module.integration.jacobiInversion  , module.coords.node, module.integration.dN);
	module.controller.addInput(module.integration.jacobiDeterminant, module.coords.node, module.integration.dN);
	module.controller.addInput(module.integration.dND              , module.coords.node, module.integration.dN);
	module.controller.prepare(module.integration.jacobiInversion, module.integration.jacobiDeterminant, module.integration.dND);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, ElementJacobian2D>(interval, module.controller, module.coords.node, module.integration.dN, module.integration.jacobiInversion, module.integration.jacobiDeterminant, module.integration.dND));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, ElementJacobian3D>(interval, module.controller, module.coords.node, module.integration.dN, module.integration.jacobiInversion, module.integration.jacobiDeterminant, module.integration.dND));
		}
	}
}

template <class Module>
void _boundaryIntegration(Module &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			module.controller.addInput(module.integration.boundary.jacobian.regions[r], module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r]);
			module.controller.prepare(module.integration.boundary.jacobian.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				if (info::mesh->boundaryRegions[r]->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryFaceJacobian>(r, interval, module.controller, module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r], module.integration.boundary.jacobian.regions[r]));
				}
				if (info::mesh->boundaryRegions[r]->dimension == 1) {
					if (info::mesh->dimension == 3) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryEdge3DJacobian>(r, interval, module.controller, module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r], module.integration.boundary.jacobian.regions[r]));
					}
					if (info::mesh->dimension == 2) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryEdge2DJacobian>(r, interval, module.controller, module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r], module.integration.boundary.jacobian.regions[r]));
					}
				}
			}
		}
	}
}

template <class Module>
void _boundaryIntegrationWithNormal(Module &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			module.controller.addInput(module.integration.boundary.jacobian.regions[r], module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r]);
			module.controller.prepare(module.integration.boundary.jacobian.regions[r]);
			module.controller.addInput(module.integration.boundary.normal.regions[r], module.integration.boundary.jacobian.regions[r]);
			module.controller.prepare(module.integration.boundary.normal.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				if (info::mesh->boundaryRegions[r]->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryFaceNormal>(r, interval, module.controller,
							module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r],
							module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r]));
				}
				if (info::mesh->boundaryRegions[r]->dimension == 1) {
//					if (info::mesh->dimension == 3) { // how to do that?
//						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryEdge3DJacobian>(r, interval, module.controller, module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r], module.integration.boundary.jacobian.regions[r]));
//					}
					if (info::mesh->dimension == 2) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, BoundaryEdge2DNormal>(r, interval, module.controller,
								module.coords.boundary.node.regions[r], module.integration.boundary.dN.regions[r],
								module.integration.boundary.jacobian.regions[r], module.integration.boundary.normal.regions[r]));
					}
				}
			}
		}
	}
}

void elementIntegration(HeatTransfer &module)
{
	_elementIntegration(module);
	_boundaryIntegration(module);
}

void elementIntegration(Acoustic &module)
{
	_elementIntegration(module);
	_boundaryIntegration(module);
}

void elementIntegration(StructuralMechanics &module)
{
	_elementIntegration(module);
	if (module.configuration.normal_pressure.size()) {
		_boundaryIntegrationWithNormal(module);
	} else {
		_boundaryIntegration(module);
	}
}

}
