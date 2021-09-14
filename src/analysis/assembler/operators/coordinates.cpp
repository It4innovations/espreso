
#include "coordinates.h"
#include "gausspoints.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "basis/expression/variable.h"
#include "config/ecf/physics/heattransfer.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template <class Module>
void _elementCoordinates(Module &module)
{
	module.controller.addInput(module.coords.node, info::mesh->nodes->coordinates);
	module.controller.prepare(module.coords.node);

	bool toGPs = Variable::list.egps.count("COORDINATE_X") || Variable::list.egps.count("COORDINATE_Y") || Variable::list.egps.count("COORDINATE_Z");
	if (toGPs) {
		module.controller.addInput(module.coords.gp, module.coords.node);
		module.controller.prepare(module.coords.gp);
	}

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(instantiate<Coordinates2DToElementNodes>(interval, module.controller, procNodes, module.coords.node));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(instantiate<Coordinates3DToElementNodes>(interval, module.controller, procNodes, module.coords.node));
		}
		if (toGPs) {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
			}
		}
	}

	auto it = Variable::list.enodes.end();
	if ((it = Variable::list.enodes.find("COORDINATE_X")) != Variable::list.enodes.end()) {
		it->second = new ParameterVariable(module.coords.node.data, module.coords.node.isconst, module.coords.node.update, 0, info::mesh->dimension);
	}
	if ((it = Variable::list.enodes.find("COORDINATE_Y")) != Variable::list.enodes.end()) {
		it->second = new ParameterVariable(module.coords.node.data, module.coords.node.isconst, module.coords.node.update, 1, info::mesh->dimension);
	}
	if (info::mesh->dimension == 3 && (it = Variable::list.enodes.find("COORDINATE_Z")) != Variable::list.enodes.end()) {
		it->second = new ParameterVariable(module.coords.node.data, module.coords.node.isconst, module.coords.node.update, 2, info::mesh->dimension);
	}

	if ((it = Variable::list.egps.find("COORDINATE_X")) != Variable::list.egps.end()) {
		it->second = new ParameterVariable(module.coords.gp.data, module.coords.gp.isconst, module.coords.gp.update, 0, info::mesh->dimension);
	}
	if ((it = Variable::list.egps.find("COORDINATE_Y")) != Variable::list.egps.end()) {
		it->second = new ParameterVariable(module.coords.gp.data, module.coords.gp.isconst, module.coords.gp.update, 1, info::mesh->dimension);
	}
	if (info::mesh->dimension == 3 && (it = Variable::list.egps.find("COORDINATE_Z")) != Variable::list.egps.end()) {
		it->second = new ParameterVariable(module.coords.gp.data, module.coords.gp.isconst, module.coords.gp.update, 2, info::mesh->dimension);
	}

	if ((it = Variable::list.node.find("COORDINATE_X")) != Variable::list.node.end()) {
		it->second = new SerializedPointsVariable(info::mesh->nodes->coordinates, 0);
	}
	if ((it = Variable::list.node.find("COORDINATE_Y")) != Variable::list.node.end()) {
		it->second = new SerializedPointsVariable(info::mesh->nodes->coordinates, 1);
	}
	if (info::mesh->dimension == 3 && (it = Variable::list.node.find("COORDINATE_Z")) != Variable::list.node.end()) {
		it->second = new SerializedPointsVariable(info::mesh->nodes->coordinates, 2);
	}
}

template <class Module>
void _boundaryCoordinates(Module &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			bool toGPs = Variable::list.region[r].egps.count("COORDINATE_X") || Variable::list.region[r].egps.count("COORDINATE_Y") || Variable::list.region[r].egps.count("COORDINATE_Z");

			module.controller.addInput(module.coords.boundary.node.regions[r], info::mesh->nodes->coordinates);
			module.controller.prepare(module.coords.boundary.node.regions[r]);

			if (toGPs) {
				module.controller.addInput(module.coords.boundary.gp.regions[r], module.coords.boundary.node.regions[r]);
				module.controller.prepare(module.coords.boundary.gp.regions[r]);
			}

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[interval].begin;
				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(instantiate<Coordinates2DToElementNodes>(interval, module.controller, procNodes, module.coords.boundary.node.regions[r]));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(instantiate<Coordinates3DToElementNodes>(interval, module.controller, procNodes, module.coords.boundary.node.regions[r]));
				}
				if (toGPs) {
					if (info::mesh->dimension == 2) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(r, interval, module.controller, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
					}
					if (info::mesh->dimension == 3) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(r, interval, module.controller, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
					}
				}
			}
		} else {
			if (Variable::list.region[r].enodes.count("COORDINATE_X") || Variable::list.region[r].enodes.count("COORDINATE_Y") || Variable::list.region[r].enodes.count("COORDINATE_Z")) {
				module.controller.addInput(module.coords.boundary.node.regions[r], info::mesh->nodes->coordinates);
				module.controller.prepare(module.coords.boundary.node.regions[r]);

				for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					auto nodes = info::mesh->boundaryRegions[r]->nodes->cbegin(t);
					if (info::mesh->dimension == 2) {
						module.boundaryOps[r][t].emplace_back(instantiate<Coordinates2DToElementNodes>(t, module.controller, nodes, module.coords.boundary.node.regions[r]));
					}
					if (info::mesh->dimension == 3) {
						module.boundaryOps[r][t].emplace_back(instantiate<Coordinates3DToElementNodes>(t, module.controller, nodes, module.coords.boundary.node.regions[r]));
					}
				}
			}
		}

		auto it = Variable::list.region[r].enodes.end();
		if ((it = Variable::list.region[r].enodes.find("COORDINATE_X")) != Variable::list.region[r].enodes.end()) {
			it->second = new ParameterVariable(module.coords.boundary.node.regions[r].data, module.coords.boundary.node.regions[r].isconst, module.coords.boundary.node.regions[r].update, 0, info::mesh->dimension);
		}
		if ((it = Variable::list.region[r].enodes.find("COORDINATE_Y")) != Variable::list.region[r].enodes.end()) {
			it->second = new ParameterVariable(module.coords.boundary.node.regions[r].data, module.coords.boundary.node.regions[r].isconst, module.coords.boundary.node.regions[r].update, 1, info::mesh->dimension);
		}
		if (info::mesh->dimension == 3 && (it = Variable::list.region[r].enodes.find("COORDINATE_Z")) != Variable::list.region[r].enodes.end()) {
			it->second = new ParameterVariable(module.coords.boundary.node.regions[r].data, module.coords.boundary.node.regions[r].isconst, module.coords.boundary.node.regions[r].update, 2, info::mesh->dimension);
		}

		if ((it = Variable::list.region[r].egps.find("COORDINATE_X")) != Variable::list.region[r].egps.end()) {
			it->second = new ParameterVariable(module.coords.boundary.gp.regions[r].data, module.coords.boundary.gp.regions[r].isconst, module.coords.boundary.gp.regions[r].update, 0, info::mesh->dimension);
		}
		if ((it = Variable::list.region[r].egps.find("COORDINATE_Y")) != Variable::list.region[r].egps.end()) {
			it->second = new ParameterVariable(module.coords.boundary.gp.regions[r].data, module.coords.boundary.gp.regions[r].isconst, module.coords.boundary.gp.regions[r].update, 1, info::mesh->dimension);
		}
		if (info::mesh->dimension == 3 && (it = Variable::list.region[r].egps.find("COORDINATE_Z")) != Variable::list.region[r].egps.end()) {
			it->second = new ParameterVariable(module.coords.boundary.gp.regions[r].data, module.coords.boundary.gp.regions[r].isconst, module.coords.boundary.gp.regions[r].update, 2, info::mesh->dimension);
		}
	}
}

void elementCoordinates(AX_HeatTransfer &module)
{
	_elementCoordinates(module);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		auto temp = module.configuration.temperature.find(info::mesh->boundaryRegions[r]->name);
		if (temp != module.configuration.temperature.end()) {
			Variable::analyze(temp->second, r);
		}

		auto flow = module.configuration.heat_flow.find(info::mesh->boundaryRegions[r]->name);
		if (flow != module.configuration.heat_flow.end()) {
			Variable::analyze(flow->second, r);
		}

		auto flux = module.configuration.heat_flux.find(info::mesh->boundaryRegions[r]->name);
		if (flux != module.configuration.heat_flux.end()) {
			Variable::analyze(flux->second, r);
		}
	}
	_boundaryCoordinates(module);
}

void elementCoordinates(AX_Acoustic &module)
{
	_elementCoordinates(module);
	_boundaryCoordinates(module);
}

}
