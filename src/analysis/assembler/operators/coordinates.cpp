
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
	module.coords.node.addInput(info::mesh->nodes->coordinates);
	module.coords.node.resize();
	module.addParameter(module.coords.node);

	bool toGPs = Variable::list.egps.count("COORDINATE_X") || Variable::list.egps.count("COORDINATE_Y") || Variable::list.egps.count("COORDINATE_Z");
	if (toGPs) {
		module.coords.gp.addInput(module.coords.node);
		module.coords.gp.resize();
		module.addParameter(module.coords.gp);
	}

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		if (info::mesh->dimension == 2) {
			module.elementOps[interval].emplace_back(new Coordinates2DToElementNodes(procNodes, module.coords.node, interval));
		}
		if (info::mesh->dimension == 3) {
			module.elementOps[interval].emplace_back(new Coordinates3DToElementNodes(procNodes, module.coords.node, interval));
		}
		if (toGPs) {
			if (info::mesh->dimension == 2) {
				module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(interval, module.integration.N, module.coords.node, module.coords.gp));
			}
			if (info::mesh->dimension == 3) {
				module.elementOps[interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(interval, module.integration.N, module.coords.node, module.coords.gp));
			}
		}
	}

	auto it = Variable::list.enodes.end();
	if ((it = Variable::list.enodes.find("COORDINATE_X")) != Variable::list.enodes.end()) {
		it->second = Variable(0, info::mesh->dimension, module.coords.node.data->datatarray().data());
	}
	if ((it = Variable::list.enodes.find("COORDINATE_Y")) != Variable::list.enodes.end()) {
		it->second = Variable(1, info::mesh->dimension, module.coords.node.data->datatarray().data());
	}
	if (info::mesh->dimension == 3 && (it = Variable::list.enodes.find("COORDINATE_Z")) != Variable::list.enodes.end()) {
		it->second = Variable(2, info::mesh->dimension, module.coords.node.data->datatarray().data());
	}

	if ((it = Variable::list.egps.find("COORDINATE_X")) != Variable::list.egps.end()) {
		it->second = Variable(0, info::mesh->dimension, module.coords.gp.data->datatarray().data());
	}
	if ((it = Variable::list.egps.find("COORDINATE_Y")) != Variable::list.egps.end()) {
		it->second = Variable(1, info::mesh->dimension, module.coords.gp.data->datatarray().data());
	}
	if (info::mesh->dimension == 3 && (it = Variable::list.egps.find("COORDINATE_Z")) != Variable::list.egps.end()) {
		it->second = Variable(2, info::mesh->dimension, module.coords.gp.data->datatarray().data());
	}
}

template <class Module>
void _boundaryCoordinates(Module &module)
{
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			bool toGPs = Variable::list.region[r].egps.count("COORDINATE_X") || Variable::list.region[r].egps.count("COORDINATE_Y") || Variable::list.region[r].egps.count("COORDINATE_Z");

			module.coords.boundary.node.regions[r].addInput(info::mesh->nodes->coordinates);
			module.coords.boundary.node.regions[r].resize();
			module.addParameter(module.coords.boundary.node.regions[r]);

			if (toGPs) {
				module.coords.boundary.gp.regions[r].addInput(module.coords.boundary.node.regions[r]);
				module.coords.boundary.gp.regions[r].resize();
				module.addParameter(module.coords.boundary.gp.regions[r]);
			}

			for(size_t interval = 0; interval < info::mesh->boundaryRegions[r]->eintervals.size(); ++interval) {
				auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[interval].begin;
				if (info::mesh->dimension == 2) {
					module.boundaryOps[r][interval].emplace_back(new Coordinates2DToElementNodes(procNodes, module.coords.boundary.node.regions[r], interval));
				}
				if (info::mesh->dimension == 3) {
					module.boundaryOps[r][interval].emplace_back(new Coordinates3DToElementNodes(procNodes, module.coords.boundary.node.regions[r], interval));
				}
				if (toGPs) {
					if (info::mesh->dimension == 2) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 2, FromNodesToGaussPoints>(r, interval, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
					}
					if (info::mesh->dimension == 3) {
						module.boundaryOps[r][interval].emplace_back(instantiate<typename Module::NGP, 3, FromNodesToGaussPoints>(r, interval, module.integration.boundary.N.regions[r], module.coords.boundary.node.regions[r], module.coords.boundary.gp.regions[r]));
					}
				}
			}

			auto it = Variable::list.region[r].enodes.end();
			if ((it = Variable::list.region[r].enodes.find("COORDINATE_X")) != Variable::list.region[r].enodes.end()) {
				it->second = Variable(0, info::mesh->dimension, module.coords.boundary.node.regions[r].data->datatarray().data());
			}
			if ((it = Variable::list.region[r].enodes.find("COORDINATE_Y")) != Variable::list.region[r].enodes.end()) {
				it->second = Variable(1, info::mesh->dimension, module.coords.boundary.node.regions[r].data->datatarray().data());
			}
			if (info::mesh->dimension == 3 && (it = Variable::list.region[r].enodes.find("COORDINATE_Z")) != Variable::list.region[r].enodes.end()) {
				it->second = Variable(2, info::mesh->dimension, module.coords.boundary.node.regions[r].data->datatarray().data());
			}

			if ((it = Variable::list.region[r].egps.find("COORDINATE_X")) != Variable::list.region[r].egps.end()) {
				it->second = Variable(0, info::mesh->dimension, module.coords.boundary.gp.regions[r].data->datatarray().data());
			}
			if ((it = Variable::list.region[r].egps.find("COORDINATE_Y")) != Variable::list.region[r].egps.end()) {
				it->second = Variable(1, info::mesh->dimension, module.coords.boundary.gp.regions[r].data->datatarray().data());
			}
			if (info::mesh->dimension == 3 && (it = Variable::list.region[r].egps.find("COORDINATE_Z")) != Variable::list.region[r].egps.end()) {
				it->second = Variable(2, info::mesh->dimension, module.coords.boundary.gp.regions[r].data->datatarray().data());
			}
		} else {
			auto it = Variable::list.node.end();
			if ((it = Variable::list.node.find("COORDINATE_X")) != Variable::list.node.end()) {
				it->second = Variable(0, info::mesh->dimension, &info::mesh->nodes->coordinates->datatarray().data()->x);
			}
			if ((it = Variable::list.node.find("COORDINATE_Y")) != Variable::list.node.end()) {
				it->second = Variable(1, info::mesh->dimension, &info::mesh->nodes->coordinates->datatarray().data()->x);
			}
			if (info::mesh->dimension == 3 && (it = Variable::list.node.find("COORDINATE_Z")) != Variable::list.node.end()) {
				it->second = Variable(2, info::mesh->dimension, &info::mesh->nodes->coordinates->datatarray().data()->x);
			}
		}
	}
}

void elementCoordinates(AX_HeatTransfer &module)
{
	_elementCoordinates(module);

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
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
