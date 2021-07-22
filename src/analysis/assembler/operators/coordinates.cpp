
#include "coordinates.h"
#include "gausspoints.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template <class Module>
void _elementCoordinates(Module &module)
{
	module.coords.node.addInput(info::mesh->nodes->coordinates);
	module.coords.gp.addInput(module.coords.node);
	module.coords.node.resize();
	module.coords.gp.resize();

	module.addParameter(module.coords.node);
	module.addParameter(module.coords.gp);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(new Coordinates2DToElementNodes(procNodes, module.coords.node, interval));
			module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(interval, module.integration.N, module.coords.node, module.coords.gp));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(new Coordinates3DToElementNodes(procNodes, module.coords.node, interval));
			module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(interval, module.integration.N, module.coords.node, module.coords.gp));
		}
	}
}

template <class Module>
void _boundaryCoordinates(Module &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			module.coords.boundary.node.regions[r].addInput(info::mesh->nodes->coordinates);
			module.coords.boundary.gp.regions[r].addInput(module.coords.boundary.node.regions[r]);
			module.coords.boundary.node.regions[r].isset = true;
			module.coords.boundary.gp.regions[r].isset = true;
			module.coords.boundary.node.regions[r].resize();
			module.coords.boundary.gp.regions[r].resize();

			module.addParameter(module.coords.boundary.node.regions[r]);
			module.addParameter(module.coords.boundary.gp.regions[r]);

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[interval].begin;
				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(new Coordinates2DToElementNodes(procNodes, module.coords.boundary.node.regions[r], interval));
					module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(r, interval, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(new Coordinates3DToElementNodes(procNodes, module.coords.boundary.node.regions[r], interval));
					module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(r, interval, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
				}
			}
		}
	}
}

void elementCoordinates(AX_HeatTransfer &module)
{
	_elementCoordinates(module);
	_boundaryCoordinates(module);
}

void elementCoordinates(AX_Acoustic &module)
{
	_elementCoordinates(module);
	_boundaryCoordinates(module);
}

//BoundaryCoordinates::BoundaryCoordinates(AX_HeatTransfer &module): BoundaryOperatorBuilder("BOUNDARY COORDINATES")
//{
//	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
//		module.coords.boundary.node.regions[r].addInput(info::mesh->nodes->coordinates);
//		module.coords.boundary.gp.regions[r].addInput(module.coords.boundary.node.regions[r]);
//		module.coords.boundary.node.regions[r].isset = true;
//		module.coords.boundary.gp.regions[r].isset = true;
//		module.coords.boundary.node.regions[r].resize();
//		module.coords.boundary.gp.regions[r].resize();
//
//		module.addParameter(module.coords.boundary.node.regions[r]);
//		module.addParameter(module.coords.boundary.gp.regions[r]);
//	}
//}
//
//void BoundaryCoordinates::apply(int region, int interval)
//{
////	auto procNodes = info::mesh->boundaryRegions[region]->nodes->cbegin() + info::mesh->boundaryRegions[region]->eintervals[interval].begin;
////	if (info::mesh->dimension == 2) {
////		iterate_boundary(Coordinates2DToElementNodes(procNodes, module.coords.boundary.node.regions[region], interval), region);
////		iterate_boundary_gps<Assembler::NGP>(FromNodesToGaussPoints<2>(module.integration.boundary.N.regions[region], module.coords.boundary.node.regions[region], module.coords.boundary.gp.regions[region], interval), region);
////	}
////	if (info::mesh->dimension == 3) {
////		iterate_boundary(Coordinates3DToElementNodes(procNodes, module.coords.boundary.node.regions[region], interval), region);
////		iterate_boundary_gps<Assembler::NGP>(FromNodesToGaussPoints<3>(module.integration.boundary.N.regions[region], module.coords.boundary.node.regions[region], module.coords.boundary.gp.regions[region], interval), region);
////	}
//}

}
