
#include "copy.h"

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/module/heattransfer.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/nodestore.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

void averageEnodesToNodes(const ParameterData &from, NodeData &to)
{
	std::fill(to.data.begin(), to.data.end(), 0);
	auto procNodes = info::mesh->elements->nodes->begin();
	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		InputParameterIterator input(from, i, to.dimension);
		for (esint e = info::mesh->elements->eintervals[i].begin; e < info::mesh->elements->eintervals[i].end; ++e, ++procNodes) {
			for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++input) {
				for (int d = 0; d < to.dimension; ++d) {
					to.data[*n * to.dimension + d] += input[d];
				}
			}
		}
	}

	auto &nelements = info::mesh->nodes->elements->boundarytarray();
	for (size_t i = 0; i < to.data.size() / to.dimension; i++) {
		for (int d = 0; d < to.dimension; d++) {
			to.data[i * to.dimension + d] /= nelements[i + 1] - nelements[i];
		}
	}

	std::vector<std::vector<double> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

	auto nranks = info::mesh->nodes->ranks->begin();
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbors[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < to.dimension; d++) {
						sBuffer[noffset].push_back(to.data[n * to.dimension + d]);
					}
				}
			}
		}
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
		rBuffer[n].resize(sBuffer[n].size());
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("exchange diagonal values.\n");
	}

	nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> nindex(info::mesh->neighbors.size());
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		if (nranks->size() > 1) {
			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				if (*r != info::mpi::rank) {
					while (info::mesh->neighbors[noffset] < *r) {
						++noffset;
					}
					for (int d = 0; d < to.dimension; d++) {
						to.data[n * to.dimension + d] += rBuffer[noffset][nindex[noffset]++];
					}
				}
			}
		}
	}
}

void copyNodesToEnodes(HeatTransfer &module, const NodeData &from, ParameterData &to)
{
	module.controller.addInput(to, info::mesh->nodes->coordinates);
	module.controller.prepare(to);
	to.setUpdate(1);

	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
		switch (from.dimension) {
		case 1: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<1> >(interval, module.controller, from, procNodes, to)); break;
		case 2: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<2> >(interval, module.controller, from, procNodes, to)); break;
		case 3: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<3> >(interval, module.controller, from, procNodes, to)); break;
		default: eslog::globalerror("ESPRESO internal error: not-implemented dimension inside: copyNodesToEnodes\n");
		}
	}
}

void copyNodesToBnodes(HeatTransfer &module, const NodeData &from, ParameterData &to, size_t region)
{
	module.controller.addInput(to, info::mesh->boundaryRegions[region]->elements);
	module.controller.prepare(to);
	to.setUpdate(1);

	for(size_t interval = 0; interval < info::mesh->boundaryRegions[region]->eintervals.size(); ++interval) {
		auto procNodes = info::mesh->boundaryRegions[region]->elements->cbegin() + info::mesh->boundaryRegions[region]->eintervals[interval].begin;
		switch (from.dimension) {
		case 1: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<1> >(interval, module.controller, from, procNodes, to)); break;
		case 2: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<2> >(interval, module.controller, from, procNodes, to)); break;
		case 3: module.elementRes[interval].emplace_back(instantiate<CopyNodesToEnodes<3> >(interval, module.controller, from, procNodes, to)); break;
		default: eslog::globalerror("ESPRESO internal error: not-implemented dimension inside: copyNodesToEnodes\n");
		}
	}
}

}
